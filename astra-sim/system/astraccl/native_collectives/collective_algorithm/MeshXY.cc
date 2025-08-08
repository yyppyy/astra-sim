/******************************************************************************
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
*******************************************************************************/

#include "astra-sim/system/astraccl/native_collectives/collective_algorithm/MeshXY.hh"

#include "astra-sim/system/PacketBundle.hh"
#include "astra-sim/system/RecvPacketEventHandlerData.hh"

using namespace AstraSim;

inline std::unordered_map<int, int> MeshXY::instance_count_;

static inline int step_id_col_major(int id,
                                 int mesh_x, int mesh_y,
                                 int delta_x, int delta_y) {
    if (id < 0 || mesh_x <= 0 || mesh_y <= 0) return -1;
    int x = id / mesh_y;   // column-major unflatten
    int y = id % mesh_y;

    int nx = x + delta_x;
    int ny = y + delta_y;
    if (nx < 0 || nx >= mesh_x || ny < 0 || ny >= mesh_y) return -1;

    return nx * mesh_y + ny;   // column-major flatten
}

static inline int get_tag_from_id(int id, int instance_id) {
    return (instance_id << 24) + id;
}

static inline int get_src_from_tag(int tag) {
    return tag % (1 << 24);
}

MeshXY::MeshXY(ComType type,
           int id,
           MeshTopology* mesh_topology,
           uint64_t data_size,
           InjectionPolicy injection_policy,
           int group_x,
           int group_y,
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

    this->group_x_ = group_x;
    this->group_y_ = group_y;
    assert(this->mesh_x_ % group_x == 0);
    assert(this->mesh_y_ % group_y == 0);
    this->group_i_ = this->mesh_i_ / group_x;
    this->group_j_ = this->mesh_i_ / group_y;

    this->part_x_ = part_x;
    this->part_y_ = part_y;
    assert(this->group_x_ % part_x == 0);
    assert(this->group_y_ % part_y == 0);
    this->part_i_ = (this->mesh_i_ % group_x) / part_x;
    this->part_j_ = (this->mesh_j_ % group_y) / part_y;

    assert(type != ComType::All_Reduce);
    this->inter_part_ = inter_part;
    // this->in_x_phase_ = true;

    this->send_ups_ = std::vector<int>{};
    this->send_downs_ = std::vector<int>{};
    this->send_lefts_ = std::vector<int>{};
    this->send_rights_ = std::vector<int>{};

    this->recv_ups_ = std::unordered_set<int>{};
    this->recv_downs_ = std::unordered_set<int>{};
    this->recv_lefts_ = std::unordered_set<int>{};
    this->recv_rights_ = std::unordered_set<int>{};

    if (inter_part) {
        int up = step_id_col_major(this->id, this->mesh_x_, this->mesh_y_, -part_x * this->part_i_, 0);
        if (up != this->id) {
            this->recv_ups_.insert(up);
        }
        int down = step_id_col_major(this->id, this->mesh_x_, this->mesh_y_, part_x * (group_x / part_x - this->part_i_ - 1), 0);
        if (down != this->id) {
            this->recv_downs_.insert(down);            
        }
        int left = step_id_col_major(this->id, this->mesh_x_, this->mesh_y_, 0, -part_y * this->part_j_);
        if (left != this->id) {
            this->recv_lefts_.insert(left);
        }
        int right = step_id_col_major(this->id, this->mesh_x_, this->mesh_y_, 0, part_y * (group_y / part_y - this->part_j_ - 1));
        if (right != this->id) {
            this->recv_rights_.insert(right);
        }
        if (this->part_i_ == 0 || this->part_i_ == ((group_x / part_x) - 1)) {
            int step = 0;
            for (int i = step_id_col_major(id, this->mesh_x_, this->mesh_y_, -part_x, 0);
                    step < this->part_i_;
                    i = step_id_col_major(i, this->mesh_x_, this->mesh_y_, -part_x, 0)) {
                this->send_ups_.push_back(i);
                ++step;
            }
            step = 0;
            for (int i = step_id_col_major(id, this->mesh_x_, this->mesh_y_, part_x, 0);
                    step < group_x / part_x - this->part_i_ - 1;
                    i = step_id_col_major(i, this->mesh_x_, this->mesh_y_, part_x, 0)) {
                // printf("debug1 id:%d, type:%d, down push:%d\n", this->id, this->comType, i);
                this->send_downs_.push_back(i);
                ++step;
            }
        }
        if (this->part_j_ == 0 || this->part_j_ == ((group_y / part_y) - 1)) {        
            int step = 0;
            for (int i = step_id_col_major(id, this->mesh_x_, this->mesh_y_, 0, -part_y);
                    step < this->part_j_;
                    i = step_id_col_major(i, this->mesh_x_, this->mesh_y_, 0, -part_y)) {
                this->send_lefts_.push_back(i);
                ++step;
            }
            step = 0;
            for (int i = step_id_col_major(id, this->mesh_x_, this->mesh_y_, 0, part_y);
                    step < group_y / part_y - this->part_j_ - 1;
                    i = step_id_col_major(i, this->mesh_x_, this->mesh_y_, 0, part_y)) {
                this->send_rights_.push_back(i);
                ++step;
            }
        }
    } else {
        int up = step_id_col_major(this->id, this->mesh_x_, this->mesh_y_, -(this->mesh_i_ % part_x), 0);
        if (up != this->id) {
            this->recv_ups_.insert(up);
        }
        int down = step_id_col_major(this->id, this->mesh_x_, this->mesh_y_, part_x - (this->mesh_i_ % part_x) - 1, 0);
        if (down != this->id) {
            this->recv_downs_.insert(down);            
        }
        int left = step_id_col_major(this->id, this->mesh_x_, this->mesh_y_, 0, -(this->mesh_j_ % part_y));
        if (left != this->id) {
            this->recv_lefts_.insert(left);
        }
        int right = step_id_col_major(this->id, this->mesh_x_, this->mesh_y_, 0, part_y - (this->mesh_j_ % part_y) - 1);
        if (right != this->id) {
            this->recv_rights_.insert(right);
        }
        if ((this->mesh_i_ % part_x == 0) || (this->mesh_i_ % part_x) == (part_x - 1)) {
            int step = 0;
            for (int i = step_id_col_major(id, this->mesh_x_, this->mesh_y_, -1, 0);
                    step < (this->mesh_i_ % part_x);
                    i = step_id_col_major(i, this->mesh_x_, this->mesh_y_, -1, 0)) {
                this->send_ups_.push_back(i);
                ++step;
            }
            step = 0;
            for (int i = step_id_col_major(id, this->mesh_x_, this->mesh_y_, 1, 0);
                    step < part_x - (this->mesh_i_ % part_x) - 1;
                    i = step_id_col_major(i, this->mesh_x_, this->mesh_y_, 1, 0)) {
                this->send_downs_.push_back(i);
                ++step;
            }
        }
        if ((this->mesh_j_ % part_y == 0) || (this->mesh_j_ % part_y) == (part_y - 1)) {
            int step = 0;
            for (int i = step_id_col_major(id, this->mesh_x_, this->mesh_y_, 0, -1);
                    step < (this->mesh_j_ % part_y);
                    i = step_id_col_major(i, this->mesh_x_, this->mesh_y_, 0, -1)) {
                this->send_lefts_.push_back(i);
                ++step;
            }
            step = 0;
            for (int i = step_id_col_major(id, this->mesh_x_, this->mesh_y_, 0, 1);
                    step < part_y - (this->mesh_j_ % part_y) - 1;
                    i = step_id_col_major(i, this->mesh_x_, this->mesh_y_, 0, 1)) {
                this->send_rights_.push_back(i);
                ++step;
            }
        }
    }

    // data size is already divided by N or K splits when fed in
    assert(data_size % (group_x * group_y) == 0);
    this->x_phase_msg_size_ = data_size / group_x;
    this->y_phase_msg_size_ = data_size / group_x * group_y;

    this->packets_to_up_sent_ = 0;
    this->packets_to_down_sent_ = 0;
    this->packets_to_left_sent_ = 0;
    this->packets_to_right_sent_ = 0;

    this->done_x_phase_send_ = false;
    this->done_x_phase_recv_ = false;
    this->done_y_phase_send_ = false;
    this->done_y_phase_recv_ = false;
    this->done_alltoall_send_ = false;
    this->done_alltoall_recv_ = false;


    /////////////////////////////////// above are for all-gather and reduce-scatter /////////////////////////////////
    ///////////////////////////////////////////// below are for all-toall ///////////////////////////////////////////

    this->alltoall_send_matrix_ = alltoall_send_matrix;
    this->alltoall_recv_matrix_ = alltoall_recv_matrix;

    this->alltoall_packet_sent_ = 0;
    this->alltoall_packet_recved_ = 0;

    ++(instance_count_[this->id]);
    this->instance_id_ = instance_count_[this->id];

    LoggerFactory::get_logger("system::collective::MeshXY")
        ->debug("id:{} instance:{}, type:{}, mesh:({},{}), group:({},{}), partition:({},{}), mesh_idx:({},{}), group_idx:({},{}), partition_idx:({},{}), inter_part:{}, x_msg_size:{}, y_msg_size:{}",
               this->id, this->instance_id_, (int)(this->comType),
                this->mesh_x_, this->mesh_y_, this->group_x_, this->group_y_, this->part_x_, this->part_y_,
                this->mesh_i_, this->mesh_j_, this->group_i_, this->group_j_, this->part_i_, this->part_j_,
                this->inter_part_, this->x_phase_msg_size_, this->y_phase_msg_size_);

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

bool MeshXY::all_done() {
    if (this->comType == ComType::All_to_All) {
        return this->done_alltoall_send_ && this->done_alltoall_recv_;
    } else {
        return this->done_x_phase_send_ && this->done_x_phase_recv_ && this->done_y_phase_send_ && this->done_y_phase_recv_;
    }
}

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
        bool send_msg = true;

        if (comType == ComType::All_to_All) {
            if (this->alltoall_send_matrix_.size() == 0) {
                // a dummy packet to wake up and exit
                send_msg = false;
            } else {
                dst = this->alltoall_send_matrix_[this->alltoall_packet_sent_].first;
                msg_size = this->alltoall_send_matrix_[this->alltoall_packet_sent_].second;
                ++(this->alltoall_packet_sent_);
            }
            if (this->alltoall_packet_sent_ == this->alltoall_send_matrix_.size()) {
                this->done_alltoall_send_ = true;
            }
        } else if (!this->done_x_phase_send_) {
            bool send2up;
            if (this->packets_to_up_sent_ == this->send_ups_.size()) {
                send2up = false;
            } else if (this->packets_to_down_sent_ == this->send_downs_.size()) {
                send2up = true;
            } else {
                send2up = this->packets_to_up_sent_ <= this->packets_to_down_sent_;
            }

            if (send2up) {
                dst = this->send_ups_[this->packets_to_up_sent_];
                ++(this->packets_to_up_sent_);
            } else {
                // printf("debug2 id: %d, type:%d, in x: %d, done x send: %d, done x recv: %d, done y send: %d, done y recv: %d, send_down idx:%d, size:%d\n", this->id, this->comType, this->in_x_phase_,
                    // this->done_x_phase_send_, this->done_x_phase_recv_, this->done_y_phase_send_, this->done_y_phase_recv_, this->packets_to_down_sent_, this->send_downs_.size());
                dst = this->send_downs_[this->packets_to_down_sent_];
                ++(this->packets_to_down_sent_);
            }
            msg_size = this->x_phase_msg_size_;

            if (this->packets_to_up_sent_ == this->send_ups_.size() && this->packets_to_down_sent_ == this->send_downs_.size()) {
                // all done with x phase send
                this->done_x_phase_send_ = true;
                // if (this->done_x_phase_recv_) {
                //     this->in_x_phase_ = false;
                // }
            }
        } else {
            bool send2left;
            if (this->packets_to_left_sent_ == this->send_lefts_.size()) {
                send2left = false;
            } else if (this->packets_to_right_sent_ == this->send_rights_.size()) {
                send2left = true;
            } else {
                send2left = this->packets_to_left_sent_ <= this->packets_to_right_sent_;
            }

            if (send2left) {
                dst = this->send_lefts_[this->packets_to_left_sent_];
                ++(this->packets_to_left_sent_);
            } else {
                // printf("%d %d\n", this->send_rights_.size(), this->packets_to_right_sent_);
                dst = this->send_rights_[this->packets_to_right_sent_];
                ++(this->packets_to_right_sent_);
            }
            msg_size = this->y_phase_msg_size_;

            if (this->packets_to_left_sent_ == this->send_lefts_.size() && this->packets_to_right_sent_ == this->send_rights_.size()) {
                // done with y phase send
                this->done_y_phase_send_ = true;
            }
        }

        if (send_msg) {
            // printf("here %d %d %ld\n", this->id, this->alltoall_packet_sent_, this->alltoall_send_matrix_.size());
            LoggerFactory::get_logger("system::collective::MeshXY")
                ->debug("id:{}, instance:{}, send packet to: {}, size: {}",
                    this->id, this->instance_id_, dst, msg_size);

            sim_request snd_req;
            // printf("src, dst %d, %d\n", this->id, dst);
            assert(stream->owner->id == this->id);
            snd_req.srcRank = this->id;
            snd_req.dstRank = dst;
            snd_req.tag = get_tag_from_id(this->id, this->instance_id_);
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
        }

        // if we don't need to recv we should exit right now
        // if ((this->comType == ComType::All_to_All && this->alltoall_recv_matrix_.size() == 0)
        //     || ((this->comType == ComType::All_Gather || this->comType == ComType::Reduce_Scatter)
        //         && !this->in_x_phase_ && this->lefts_.size() == 0 && this->rights_.size() == 0)) {
        // if (this->comType == ComType::All_to_All && this->alltoall_recv_matrix_.size() == 0) {
        if (this->all_done()) {
            LoggerFactory::get_logger("system::collective::MeshXY")
                ->debug("id:{}, skip recv and exit", this->id);
            exit();
        }

    } else if (event == EventType::PacketReceived) {

        auto ehd = (RecvPacketEventHandlerData*)data;
        int tag = ehd->tag;
        int src = get_src_from_tag(tag);

        // handle recved packets
        if (comType == ComType::All_to_All) {
            ++(this->alltoall_packet_recved_);
            LoggerFactory::get_logger("system::collective::MeshXY")
                ->debug("id:{}, instance:{}, All_to_All recv packet from: {}, total recved: {}, exp recved: {}",
                    this->id, this->instance_id_, src, this->alltoall_packet_recved_, this->alltoall_recv_matrix_.size());
            if (this->alltoall_packet_recved_ == this->alltoall_recv_matrix_.size()) {
                LoggerFactory::get_logger("system::collective::MeshXY")
                     ->debug("id:{}, done all_to_all recv and exit", this->id);
                this->done_alltoall_recv_ = true;
            }
        } else if (!this->done_x_phase_recv_) {
            bool from_up = this->recv_ups_.find(src) != this->recv_ups_.end();
            bool from_down = this->recv_downs_.find(src) != this->recv_downs_.end();
            if (from_up) {
                assert(!from_down);
                this->recv_ups_.erase(src);
            } else if (from_down) {
                assert(!from_up);
                this->recv_downs_.erase(src);
            } else{
                printf("id:%d instance:%d unexpected recv from %d\n", this->id, this->instance_id_, src);
                assert(0);
            }
            LoggerFactory::get_logger("system::collective::MeshXY")
                ->debug("id:{}, instance:{}, X-phase recv packet from: {}({}), dir to recv: {}, exp dir recved: {}",
                    id, this->instance_id_, src, from_up ? "up" : "down",
                    from_up ? this->recv_ups_.size() : this->recv_downs_.size(),
                    1);
            if (this->recv_ups_.size() == 0 && this->recv_downs_.size() == 0) {
                this->done_x_phase_recv_ = true;
                // if (this->done_x_phase_send_) {
                //     this->in_x_phase_ = false;
                // }
                
                // if (this->recv_lefts_.size() == 0 && this->recv_rights_.size() == 0) {
                //     LoggerFactory::get_logger("system::collective::MeshXY")
                //         ->debug("id:{}, skip Y-phase and exit", this->id);
                //     exit();
                //     return;
                // }

                // post recv for y phase
                std::vector<int> recv_srcs;
                for (const int left: this->recv_lefts_) {
                    recv_srcs.push_back(left);
                }
                for (const int right: this->recv_rights_) {
                    recv_srcs.push_back(right);
                }
                if (recv_srcs.size() == 0) {
                    this->done_y_phase_send_ = true;
                    this->done_y_phase_recv_ = true;
                }

                LoggerFactory::get_logger("system::collective::MeshXY")
                    ->debug("id:{}, instance:{}, post recv Y-phase packets, len: {}", this->id, this->instance_id_, recv_srcs.size());
                for (int i = 0; i < recv_srcs.size(); ++i) {
                    int recv_src = recv_srcs[i];

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
                        this->y_phase_msg_size_,
                        UINT8,
                        recv_src,
                        get_tag_from_id(recv_src, this->instance_id_),
                        &rcv_req,
                        Sys::FrontEndSendRecvType::COLLECTIVE,
                        &Sys::handleEvent,
                        ehd
                    );
                }

                // insert initial packets for the Y phase
                for (int i = 0; i < this->send_lefts_.size() + this->send_rights_.size(); ++i) {
                    (new PacketBundle(
                        stream->owner,
                        stream,
                        false, // no need to process, initial packets
                        false, // always set send_back to false for now
                        this->y_phase_msg_size_,
                        this->transmition_
                    ))->send_to_MA();
                }
                if (this->send_lefts_.size() + this->send_rights_.size() == 0) {
                    // no packets to send
                    this->done_y_phase_send_ = true;
                }
                LoggerFactory::get_logger("system::collective::MeshXY")
                    ->debug("id:{}, insert Y-phase init packets, len: {}", this->id, this->send_lefts_.size() + this->send_rights_.size());
            }
        } else {
            bool from_left = this->recv_lefts_.find(src) != this->recv_lefts_.end();
            bool from_right = this->recv_rights_.find(src) != this->recv_rights_.end();
            if (from_left) {
                assert(!from_right);
                this->recv_lefts_.erase(src);
            } else if (from_right) {
                assert(!from_left);
                this->recv_rights_.erase(src);
            } else{
                assert(0);
            }
            LoggerFactory::get_logger("system::collective::MeshXY")
                ->debug("id:{}, instance:{}, Y-phase recv packet from: {}({}), dir to recv: {}, exp dir recved: {}",
                    id, this->instance_id_, src, from_left ? "left" : "right",
                    from_left ? this->recv_lefts_.size() : this->recv_rights_.size(),
                    1);
            if (this->recv_lefts_.size() == 0 && this->recv_rights_.size() == 0) {
                LoggerFactory::get_logger("system::collective::MeshXY")
                    ->debug("id:{}, done Y-phase", this->id);
                // exit();
                this->done_y_phase_recv_ = true;
            }
        }

        // printf("debug5 id: %d, instance:%d, in x: %d, done x send: %d, done x recv: %d, done y send: %d, done y recv: %d\n", this->id, this->instance_id_, this->in_x_phase_,
                // this->done_x_phase_send_, this->done_x_phase_recv_, this->done_y_phase_send_, this->done_y_phase_recv_);

        if (this->all_done()) {
            exit();
        }

        // forward packet to NPU. fix me this trigger bug
        // int msg_size;
        // if (comType == ComType::All_to_All || !this->in_x_phase_) {
        //     // todo: fix dummy msg size for alltoall
        //     msg_size = this->y_phase_msg_size_;
        // } else {
        //     msg_size = this->x_phase_msg_size_;
        // }
        // (new PacketBundle(stream->owner, stream, false, false, msg_size,
        //                   MemBus::Transmition::Usual))
        //     ->send_to_NPU();

        // called when recved a packet
        // total_packets_received++;
        // insert_packet(nullptr);
        // todo: inject packets based on protocol state machine of X-Y reduce-scatter or all-gather
    } else if (event == EventType::StreamInit) {
        // called upon init

        // post recv
        std::vector<int> recv_srcs;
        std::vector<int> msg_sizes;

        if (comType == ComType::All_to_All) {
            for (int i = 0; i < this->alltoall_recv_matrix_.size(); ++i) {
                recv_srcs.push_back(this->alltoall_recv_matrix_[i].first);
                msg_sizes.push_back(this->alltoall_recv_matrix_[i].second);
            }
            if (this->alltoall_recv_matrix_.size() == 0) {
                this->done_alltoall_recv_ = true;
            }
        } else {
            // assert(this->in_x_phase_);
            auto& srcs1 = this->recv_ups_;
            auto& srcs2 = this->recv_downs_;
            int msg_size = this->x_phase_msg_size_;
            // printf("debug3 id:%d up len %d, down len %d\n", this->id, srcs1.size(), srcs2.size());
            if (srcs1.size() + srcs2.size() == 0) {
                // for any tile, if no packets to recv in X phase,
                // the only possibility is that the x dim is 1
                // so we can fast forward to Y phase
                LoggerFactory::get_logger("system::collective::MeshXY")
                    ->debug("id:{}, fast-forward to Y-phase", this->id);
                // this->in_x_phase_ = false;
                srcs1 = this->recv_lefts_;
                srcs2 = this->recv_rights_;
                assert(srcs1.size() + srcs2.size() > 0);
                msg_size = this->y_phase_msg_size_;
                this->done_x_phase_send_ = true;
                this->done_x_phase_recv_ = true;
            }
            for (const int src1: srcs1) {
                recv_srcs.push_back(src1);
                msg_sizes.push_back(msg_size);
            }
            for (const int src2: srcs2) {
                recv_srcs.push_back(src2);
                msg_sizes.push_back(msg_size);
            }
        }
        LoggerFactory::get_logger("system::collective::MeshXY")
            ->debug("id:{}, post recv init packets, len: {}", this->id, recv_srcs.size());
        for (int i = 0; i < recv_srcs.size(); ++i) {
            int recv_src = recv_srcs[i];
            int msg_size = msg_sizes[i];

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
                msg_size,
                UINT8,
                recv_src,
                get_tag_from_id(recv_src, this->instance_id_),
                &rcv_req,
                Sys::FrontEndSendRecvType::COLLECTIVE,
                &Sys::handleEvent,
                ehd
            );
        }

        // insert intial packets
        msg_sizes.clear();
        if (comType == ComType::All_to_All) {
            for (int i = 0; i < this->alltoall_send_matrix_.size(); ++i) {
                msg_sizes.push_back(this->alltoall_send_matrix_[i].second);
            }
            if (this->alltoall_send_matrix_.size() == 0) {
                // inject a dummy packet so that we will get called and exit
                // printf("id %d insert dummy packet\n", this->id);
                LoggerFactory::get_logger("system::collective::MeshXY")
                    ->debug("id:{}, inject dummy packet to wakeup and exit", this->id);
                msg_sizes.push_back(8);
            }
        } else {
            if (!this->done_x_phase_send_) {
            // if (this->in_x_phase_) {
                for (int i = 0; i < this->send_ups_.size() + this->send_downs_.size(); ++i) {
                    msg_sizes.push_back(this->x_phase_msg_size_);
                }
                if (msg_sizes.size() == 0) {
                    // no packets to send
                    this->done_x_phase_send_ = true;
                }
            } else {
                for (int i = 0; i < this->send_lefts_.size() + this->send_rights_.size(); ++i) {
                    msg_sizes.push_back(this->y_phase_msg_size_);
                }           
            }
        }
        LoggerFactory::get_logger("system::collective::MeshXY")
            ->debug("id:{}, insert init packets, len: {}", this->id, msg_sizes.size());
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
