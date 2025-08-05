/******************************************************************************
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
*******************************************************************************/

#include "astra-sim/system/astraccl/native_collectives/collective_algorithm/MeshXY.hh"

#include "astra-sim/system/PacketBundle.hh"
#include "astra-sim/system/RecvPacketEventHandlerData.hh"

using namespace AstraSim;

MeshXY::MeshXY(ComType type,
           int id,
           MeshTopology* mesh_topology,
           uint64_t data_size,
           InjectionPolicy injection_policy,
           int part_x,
           int part_y,
           bool inter_part,
           std::vector<std::pair<int, int>> alltoall_send_matrix,
           std::vector<std::pair<int, int>> alltoall_recv_matrix)
    : Algorithm() {

    this->name = Name::MeshXY;
    this->comType = type;
    this->id = id;
    this->logical_topo = mesh_topology;

    this->mesh_x_ = mesh_topology->get_x();
    this->mesh_y_ = mesh_topology->get_y();
    this->mesh_i_ = id / this->mesh_y_;
    this->mesh_j_ = id % this->mesh_y_;

    this->part_x_ = part_x;
    this->part_y_ = part_y;
    assert(this->mesh_x_ % part_x == 0);
    assert(this->mesh_y_ % part_y == 0);
    this->part_i_ = id / (this->mesh_x_ / part_x);
    this->part_j_ = id / (this->mesh_y_ / part_y);

    assert(type != ComType::All_Reduce);
    this->inter_part_ = inter_part;
    this->in_x_phase_ = true;

    this->ups_ = std::vector<int>{};
    this->downs_ = std::vector<int>{};
    this->lefts_ = std::vector<int>{};
    this->rights_ = std::vector<int>{};

    this->ups_set_ = std::unordered_set<int>{};
    this->downs_set_ = std::unordered_set<int>{};
    this->lefts_set_ = std::unordered_set<int>{};
    this->rights_set_ = std::unordered_set<int>{};

    if (inter_part) {
        int step = 0;
        for (int i = id - this->mesh_y_ * part_x; step < this->part_i_; i -= this->mesh_y_ * part_x) {
            this->ups_.push_back(i);
            this->ups_set_.insert(i);
            ++step;
        }
        step = 0;
        for (int i = id + this->mesh_y_ * part_x; step < this->mesh_x_ / part_x - this->part_i_ - 1; i += this->mesh_y_ * part_x) {
            this->downs_.push_back(i);
            this->downs_set_.insert(i);
            ++step;
        }
        step = 0;
        for (int i = id - part_y; step < this->part_j_; i -= part_y) {
            this->lefts_.push_back(i);
            this->lefts_set_.insert(i);
            ++step;
        }
        step = 0;
        for (int i = id + part_y; step < this->mesh_y_ / part_y - this->part_j_ - 1; i += part_y) {
            this->rights_.push_back(i);
            this->rights_set_.insert(i);
            ++step;
        }
        // this->up = (this->part_i_ == 0) ? (-1) : (id - this->mesh_y_ * part_x);
        // this->down = (this->part_i_ == (this->mesh_x_ / part_x - 1)) ? (-1) : (id + this->mesh_y_ * part_x);
        // this->left = (this->part_j_ == 0) ? (-1) : (id - part_y);
        // this->right = (this->part_j_ == (this->mesh_y_ / part_y - 1)) ? (-1) : (id + part_y); 
    } else {
        int step = 0;
        for (int i = id - this->mesh_y_; step < (this->mesh_i_ % part_x); i -= this->mesh_y_) {
            this->ups_.push_back(i);
            this->ups_set_.insert(i);
            ++step;
        }
        step = 0;
        for (int i = id + this->mesh_y_; step < part_x - (this->mesh_i_ % part_x) - 1; i += this->mesh_y_) {
            this->downs_.push_back(i);
            this->downs_set_.insert(i);
            ++step;
        }
        step = 0;
        for (int i = id - 1; step < (this->mesh_j_ % part_y); --i) {
            this->lefts_.push_back(i);
            this->lefts_set_.insert(i);
            ++step;
        }
        step = 0;
        for (int i = id + 1; step < part_y - (this->mesh_j_ % part_y) - 1; ++i) {
            this->rights_.push_back(i);
            this->rights_set_.insert(i);
            ++step;
        }
        // this->up = (this->mesh_i_ % part_x == 0) ? (-1) : (id - this->mesh_y_);
        // this->down = (this->mesh_i_ % part_x == (part_x - 1)) ? (-1) : (id + this->mesh_y_);
        // this->left = (this->mesh_j_ % part_y == 0) ? (-1) : (id - 1);
        // this->right = (this->mesh_j_ % part_y == (part_y - 1)) ? (-1) : (id + 1);
    }

    // data size is already divided by N or K splits when fed in
    assert(data_size % (this->mesh_x_ * this->mesh_y_) == 0);
    this->x_phase_msg_size_ = data_size / this->mesh_x_;
    this->y_phase_msg_size_ = data_size / (this->mesh_x_ * this->mesh_y_);

    this->packets_to_up_sent_ = 0;
    this->packets_to_down_sent_ = 0;
    this->packets_to_left_sent_ = 0;
    this->packets_to_right_sent_ = 0;

    this->alltoall_send_matrix_ = alltoall_send_matrix;
    this->alltoall_recv_matrix_ = alltoall_recv_matrix;

    this->alltoall_packet_sent_ = 0;
    this->alltoall_packet_recved_ = 0;

    // if (inter_part) {
    //     this->x_phase_exp_packets_up2down_ = this->part_i_;
    //     this->x_phase_exp_packets_down2up_ = part_x - this->part_i_ - 1;
    //     this->y_phase_exp_packets_left2right_ = this->part_j_;
    //     this->y_phase_exp_packets_right2left_ = part_y - this->part_j_ - 1;
    // } else {
    //     this->x_phase_exp_packets_up2down_ = this->mesh_i_ % part_x;
    //     this->x_phase_exp_packets_down2up_ = part_x - (this->mesh_i_ % part_x) - 1;
    //     this->y_phase_exp_packets_left2right_ = this->mesh_j_ % part_y;
    //     this->y_phase_exp_packets_right2left_ = part_y - (this->mesh_j_ % part_y) - 1;       
    // }
    // this->x_phase_recved_packets_up2down_ = 0;
    // this->x_phase_recved_packets_down2up_ = 0;
    // this->y_phase_recved_packets_left2right_ = 0;
    // this->y_phase_recved_packets_right2left_ = 0;

    // this->x_phase_init_packets_sent_2up_ = 0;
    // this->x_phase_init_packets_sent_2down_ = 0;
    // this->y_phase_init_packets_sent_2left = 0;
    // this->y_phase_init_packets_sent_2right = 0;
    // this->nodes_in_ring = ring_topology->get_nodes_in_ring();
    // this->curr_receiver = ring_topology->get_receiver(id, direction);
    // this->curr_sender = ring_topology->get_sender(id, direction);
    // this->parallel_reduce = 1;
    this->injection_policy_ = injection_policy;
    // this->total_packets_sent = 0;
    // this->total_packets_received = 0;
    // this->free_packets = 0;
    // this->zero_latency_packets = 0;
    // this->non_zero_latency_packets = 0;
    // this->toggle = false;
    // if (mesh_topology->get_dimension() == MeshTopology::Dimension::Local) {
    //     transmition = MemBus::Transmition::Fast;
    // } else {
    this->transmition_ = MemBus::Transmition::Usual;
    // }
    // switch (type) {
    // case ComType::All_Reduce:
    //     assert(0);
    //     break;
    // case ComType::All_to_All:
    //     // todo: set stream count.
    //     // stream seems to be per chunk, per system (npu)
    //     // this entire algo also seems to be per chunk, per system

    //     // this->stream_count = ((nodes_in_ring - 1) * nodes_in_ring) / 2;
    //     // switch (injection_policy) {
    //     // case InjectionPolicy::Aggressive:
    //     //     this->parallel_reduce = nodes_in_ring - 1;
    //     //     break;
    //     // case InjectionPolicy::Normal:
    //     //     this->parallel_reduce = 1;
    //     //     break;
    //     // default:
    //     //     this->parallel_reduce = 1;
    //     //     break;
    //     // }
    //     break;
    // default:
    //     // todo: stream count seems to be the rounds of reduce the algo needs to do...
    //     // each time it does a send/recv, it basically decrease the stream count...
    //     // but wouldn't this wrong for alltoall?
    //     // stream_count = nodes_in_ring - 1;
    // }
    // if (comType ==ComType::All_to_All || comType ==ComType::All_Gather) {
    //     max_count = 0;
    // } else {
    //     max_count = nodes_in_ring - 1;
    // }
    // remained_packets_per_message = 1;
    // remained_packets_per_max_count = 1;
    
    // // todo: set msg size
    // switch (type) {
    // case ComType::All_Reduce:
    //     // this->final_data_size = data_size;
    //     // this->msg_size = data_size / nodes_in_ring;
    //     assert(0);
    //     break;
    // case ComType::All_Gather:
    //     // this->final_data_size = data_size * nodes_in_ring;
    //     // this->msg_size = data_size;
    //     break;
    // case ComType::Reduce_Scatter:
    //     // this->final_data_size = data_size / nodes_in_ring;
    //     // this->msg_size = data_size / nodes_in_ring;
    //     break;
    // case ComType::All_to_All:
    //     // this->final_data_size = data_size;
    //     // this->msg_size = data_size / nodes_in_ring;
    //     break;
    // default:;
    // }
}

// int Ring::get_non_zero_latency_packets() {
//     return (nodes_in_ring - 1) * parallel_reduce * 1;
// }

void MeshXY::run(EventType event, CallData* data) {
    // this is the only signature exported to above stack
    // the flow of sending a packet is call Send_to_MA -> trigger a eventype::General -> call front_end_sim_send
    // the flow of reciving a packet is call front_send_sim_recv -> trigger a eventype: PacketReceived -> call Send_to_NPU

    if (event == EventType::General) {
        // seems to be called after a packet is sent from NPU to MA with a end-point delay, likely representing the packet leaves the NPU?
        // free_packets += 1;
        // ready actually send out the packet, register a recv, and do a local reduce
        // ready();
        // iterable exit the algo if the algo finish
        // iteratable();
        int dst;
        int msg_size;

        if (comType == ComType::All_to_All) {
            dst = this->alltoall_send_matrix_[this->alltoall_packet_sent_].first;
            msg_size = this->alltoall_send_matrix_[this->alltoall_packet_sent_].second;
            ++(this->alltoall_packet_sent_);
        } else if (this->in_x_phase_) {
            bool send2up;
            if (this->packets_to_up_sent_ == this->ups_.size()) {
                send2up = false;
            } else if (this->packets_to_down_sent_ == this->downs_.size()) {
                send2up = true;
            } else {
                send2up = this->packets_to_up_sent_ <= this->packets_to_down_sent_;
            }

            if (send2up) {
                dst = this->ups_[this->packets_to_up_sent_];
                ++(this->packets_to_up_sent_);
            } else {
                dst = this->downs_[this->packets_to_down_sent_];
                ++(this->packets_to_down_sent_);
            }
            msg_size = this->x_phase_msg_size_;
        } else {
            bool send2left;
            if (this->packets_to_left_sent_ == this->lefts_.size()) {
                send2left = false;
            } else if (this->packets_to_right_sent_ == this->rights_.size()) {
                send2left = true;
            } else {
                send2left = this->packets_to_left_sent_ <= this->packets_to_right_sent_;
            }

            if (send2left) {
                dst = this->lefts_[this->packets_to_left_sent_];
                ++(this->packets_to_left_sent_);
            } else {
                dst = this->rights_[this->packets_to_right_sent_];
                ++(this->packets_to_right_sent_);
            }
            msg_size = this->y_phase_msg_size_;
        }

        sim_request snd_req;
        assert(stream->owner->id == this->id);
        snd_req.srcRank = this->id;
        snd_req.dstRank = dst;
        snd_req.tag = this->id;
        snd_req.reqType = UINT8;
        snd_req.vnet = this->stream->current_queue_id;
        stream->owner->front_end_sim_send(
            0,
            Sys::dummy_data,
            msg_size,
            UINT8,
            snd_req.dstRank,
            snd_req.tag,
            &snd_req,
            Sys::FrontEndSendRecvType::COLLECTIVE,
            &Sys::handleEvent,
            nullptr
        );

    } else if (event == EventType::PacketReceived) {
        
        // forward packet to NPU
        int msg_size;
        if (comType == ComType::All_to_All || !this->in_x_phase_) {
            // todo: fix dummy msg size
            msg_size = this->y_phase_msg_size_;
        } else {
            msg_size = this->x_phase_msg_size_;
        }
        (new PacketBundle(stream->owner, stream, false, false, msg_size,
                          MemBus::Transmition::Usual))
            ->send_to_NPU();
        
        
        // handle recved packets
        if (comType == ComType::All_to_All) {
            ++(this->alltoall_packet_recved_);
            if (this->alltoall_packet_recved_ == this->alltoall_recv_matrix_.size()) {
                exit();
            }
        } else if (this->in_x_phase_) {
            auto ehd = (RecvPacketEventHandlerData*)data;
            int tag = ehd->tag;
            bool from_up = this->ups_set_.find(tag) != this->ups_set_.end();
            bool from_down = this->downs_set_.find(tag) != this->downs_set_.end();
            if (from_up) {
                assert(!from_down);
                this->ups_set_.erase(tag);
            } else if (from_down) {
                assert(!from_up);
                this->downs_set_.erase(tag);
            } else{
                assert(0);
            }
            if (this->ups_set_.size() == 0 && this->downs_set_.size() == 0) {
                this->in_x_phase_ = false;
                
                // insert initial packets for the Y phase
                for (const auto& left_id : this->lefts_) {
                    (new PacketBundle(
                        stream->owner,
                        stream,
                        false, // no need to process, initial packets
                        false, // always set send_back to false for now
                        this->y_phase_msg_size_,
                        this->transmition_
                    ))->send_to_MA();
                }

                for (const auto& right_id : this->rights_) {
                    (new PacketBundle(
                        stream->owner,
                        stream,
                        false, // no need to process, initial packets
                        false, // always set send_back to false for now
                        this->y_phase_msg_size_,
                        this->transmition_
                    ))->send_to_MA();
                }
            }
        } else {
            auto ehd = (RecvPacketEventHandlerData*)data;
            int tag = ehd->tag;
            bool from_left = this->lefts_set_.find(tag) != this->lefts_set_.end();
            bool from_right = this->rights_set_.find(tag) != this->rights_set_.end();
            if (from_left) {
                assert(!from_right);
                this->lefts_set_.erase(tag);
            } else if (from_right) {
                assert(!from_left);
                this->rights_set_.erase(tag);
            } else{
                assert(0);
            }
            if (this->lefts_set_.size() == 0 && this->rights_set_.size() == 0) {
                exit();
            }
        }
        // called when recved a packet
        // total_packets_received++;
        // insert_packet(nullptr);
        // todo: inject packets based on protocol state machine of X-Y reduce-scatter or all-gather
    } else if (event == EventType::StreamInit) {
        // called upon init

        // post recv
        std::vector<int> recv_srcs;
        if (comType == ComType::All_to_All) {
            for (int i = 0; i < this->alltoall_recv_matrix_.size(); ++i) {
                recv_srcs.push_back(this->alltoall_recv_matrix_[i].first);
            }
        } else {
            assert(this->in_x_phase_);
            for (int i = 0; i < this->ups_.size(); ++i) {
                recv_srcs.push_back(this->ups_[i]);
            }
            for (int i = 0; i < this->downs_.size(); ++i) {
                recv_srcs.push_back(this->downs_[i]);
            }
        }
        for (const auto& recv_src : recv_srcs) {
            sim_request rcv_req;
            rcv_req.vnet = this->stream->current_queue_id;
            RecvPacketEventHandlerData* ehd = new RecvPacketEventHandlerData(
                stream,
                stream->owner->id,
                EventType::PacketReceived,
                stream->current_queue_id,
                stream->stream_id,
                recv_src
            );
            stream->owner->front_end_sim_recv(
                0,
                Sys::dummy_data,
                data_size,
                UINT8,
                recv_src,
                recv_src,
                &rcv_req,
                Sys::FrontEndSendRecvType::COLLECTIVE,
                &Sys::handleEvent,
                ehd
            );
        }

        // insert intial packets
        std::vector<int> msg_sizes;
        if (comType == ComType::All_to_All) {
            for (int i = 0; i < this->alltoall_send_matrix_.size(); ++i) {
                msg_sizes.push_back(this->alltoall_send_matrix_[i].second);
            }
        } else {
            assert(this->in_x_phase_);
            for (int i = 0; i < this->ups_.size() + this->downs_.size(); ++i) {
                msg_sizes.push_back(this->x_phase_msg_size_);
            }
        }
        for (const auto& msg_size : msg_sizes) {
            (new PacketBundle(
                stream->owner,
                stream,
                false, // no need to process, initial packets
                false, // always set send_back to false for now
                msg_size,
                this->transmition_
            ))->send_to_MA();
        }

    }
}

// void Ring::release_packets() {
//     for (auto packet : locked_packets) {
//         packet->set_notifier(this);
//     }
//     if (NPU_to_MA == true) {
//         (new PacketBundle(stream->owner, stream, locked_packets, processed,
//                           send_back, msg_size, transmition))
//             ->send_to_MA();
//     } else {
//         (new PacketBundle(stream->owner, stream, locked_packets, processed,
//                           send_back, msg_size, transmition))
//             ->send_to_NPU();
//     }
//     locked_packets.clear();
// }

// void Ring::process_stream_count() {
//     if (remained_packets_per_message > 0) {
//         remained_packets_per_message--;
//     }
//     if (id == 0) {
//     }
//     if (remained_packets_per_message == 0 && stream_count > 0) {
//         stream_count--;
//         if (stream_count > 0) {
//             remained_packets_per_message = 1;
//         }
//     }
//     if (remained_packets_per_message == 0 && stream_count == 0 &&
//         stream->state != StreamState::Dead) {
//         stream->changeState(StreamState::Zombie);
//     }
// }

// void Ring::process_max_count() {
//     if (remained_packets_per_max_count > 0) {
//         remained_packets_per_max_count--;
//     }
//     if (remained_packets_per_max_count == 0) {
//         max_count--;
//         release_packets();
//         remained_packets_per_max_count = 1;
//     }
// }

// void Ring::reduce() {
//     process_stream_count();
//     packets.pop_front();
//     free_packets--;
//     total_packets_sent++;
// }

// bool Ring::iteratable() {
//     if (stream_count == 0 &&
//         free_packets == (parallel_reduce * 1)) {  // && not_delivered==0
//         exit();
//         return false;
//     }
//     return true;
// }

// void Ring::insert_packet(Callable* sender) {
//     // 0-latency packets are those who has no latency for processing...meaning that the packet is initiated by self (i.g., the first packets in each parallel reduce ring)
//     // non-0-latency packets are those whose has end-point delay...meaning that the packet is forwarded by others
//     if (zero_latency_packets == 0 && non_zero_latency_packets == 0) {
//         zero_latency_packets = parallel_reduce * 1;
//         non_zero_latency_packets =
//             get_non_zero_latency_packets();  //(nodes_in_ring-1)*parallel_reduce*1;
//         toggle = !toggle;
//     }
//     if (zero_latency_packets > 0) {
//         packets.push_back(MyPacket(
//             stream->current_queue_id, curr_sender,
//             curr_receiver));  // vnet Must be changed for alltoall topology
//         packets.back().sender = sender;
//         locked_packets.push_back(&packets.back());
//         processed = false;
//         send_back = false;
//         NPU_to_MA = true;
//         process_max_count();
//         zero_latency_packets--;
//         return;
//     } else if (non_zero_latency_packets > 0) {
//         packets.push_back(MyPacket(
//             stream->current_queue_id, curr_sender,
//             curr_receiver));  // vnet Must be changed for alltoall topology
//         packets.back().sender = sender;
//         locked_packets.push_back(&packets.back());
//         if (comType == ComType::Reduce_Scatter ||
//             (comType == ComType::All_Reduce && toggle)) {
//             processed = true;
//         } else {
//             processed = false;
//         }
//         if (non_zero_latency_packets <= parallel_reduce * 1) {
//             send_back = false;
//         } else {
//             send_back = true;
//         }
//         NPU_to_MA = false;
//         process_max_count();
//         non_zero_latency_packets--;
//         return;
//     }
//     Sys::sys_panic("should not inject nothing!");
// }

// bool Ring::ready() {
//     if (stream->state == StreamState::Created ||
//         stream->state == StreamState::Ready) {
//         stream->changeState(StreamState::Executing);
//     }
//     if (packets.size() == 0 || stream_count == 0 || free_packets == 0) {
//         return false;
//     }
//     MyPacket packet = packets.front();
//     sim_request snd_req;
//     snd_req.srcRank = id;
//     snd_req.dstRank = packet.preferred_dest;
//     snd_req.tag = stream->stream_id;
//     snd_req.reqType = UINT8;
//     snd_req.vnet = this->stream->current_queue_id;
//     stream->owner->front_end_sim_send(
//         0, Sys::dummy_data, msg_size, UINT8, packet.preferred_dest,
//         stream->stream_id, &snd_req, Sys::FrontEndSendRecvType::COLLECTIVE,
//         &Sys::handleEvent,
//         nullptr);  // stream_id+(packet.preferred_dest*50)
//     sim_request rcv_req;
//     rcv_req.vnet = this->stream->current_queue_id;
//     RecvPacketEventHandlerData* ehd = new RecvPacketEventHandlerData(
//         stream, stream->owner->id, EventType::PacketReceived,
//         packet.preferred_vnet, packet.stream_id);
//     stream->owner->front_end_sim_recv(
//         0, Sys::dummy_data, msg_size, UINT8, packet.preferred_src,
//         stream->stream_id, &rcv_req, Sys::FrontEndSendRecvType::COLLECTIVE,
//         &Sys::handleEvent,
//         ehd);  // stream_id+(owner->id*50)
//     reduce();
//     return true;
// }

// void Ring::exit() {
//     if (packets.size() != 0) {
//         packets.clear();
//     }
//     if (locked_packets.size() != 0) {
//         locked_packets.clear();
//     }
//     stream->owner->proceed_to_next_vnet_baseline((StreamBaseline*)stream);
//     return;
// }

void MeshXY::exit() {
    stream->owner->proceed_to_next_vnet_baseline((StreamBaseline*)stream);
}
