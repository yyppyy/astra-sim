/******************************************************************************
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
*******************************************************************************/

#ifndef __MESHXY_HH__
#define __MESHXY_HH__

#include "astra-sim/system/MemBus.hh"
#include "astra-sim/system/MyPacket.hh"
#include "astra-sim/system/astraccl/Algorithm.hh"
#include "astra-sim/system/astraccl/native_collectives/logical_topology/MeshTopology.hh"
#include "astra-sim/common/Logging.hh"

#include <unordered_set>
#include <vector>

namespace AstraSim {

class MeshXY : public Algorithm {
  public:
    MeshXY(ComType type,
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
         std::vector<std::pair<int, int>> alltoall_recv_matrix);

    virtual void run(EventType event, CallData* data);

    // void process_stream_count();
    // void release_packets();
    // virtual void process_max_count();
    // void reduce();
    // bool iteratable();
    // virtual int get_non_zero_latency_packets();
    // void insert_packet(Callable* sender);
    // bool ready();
    void exit();

    // RingTopology::Direction dimension;
    // RingTopology::Direction direction;
    MemBus::Transmition transmition_;
    // int zero_latency_packets;
    // int non_zero_latency_packets;
    // int id;
    // int curr_receiver;
    // int curr_sender;
    // int nodes_in_ring;
    // int stream_count;
    // int max_count;
    // int remained_packets_per_max_count;
    // int remained_packets_per_message;
    // int parallel_reduce;
    InjectionPolicy injection_policy_;
    // std::list<MyPacket> packets;
    // bool toggle;
    // long free_packets;
    // long total_packets_sent;
    // long total_packets_received;
    // uint64_t msg_size;
    // std::list<MyPacket*> locked_packets;
    // bool processed;
    // bool send_back;
    // bool NPU_to_MA;

    int mesh_x_; // mesh dim
    int mesh_y_;
    int mesh_i_; // global index in the mesh
    int mesh_j_;

    int group_x_; // EP group x dim
    int group_y_;
    int group_i_; // index of group in mesh
    int group_j_;

    int part_x_; // partition x dim
    int part_y_;
    int part_i_; // index of partition in group
    int part_j_;

    bool inter_part_;
    
    bool in_x_phase_;

    std::vector<int> ups_;
    std::vector<int> downs_;
    std::vector<int> lefts_;
    std::vector<int> rights_;
    
    std::unordered_set<int> ups_set_;
    std::unordered_set<int> downs_set_;
    std::unordered_set<int> lefts_set_;
    std::unordered_set<int> rights_set_;
    
    int x_phase_msg_size_;
    int y_phase_msg_size_;

    int packets_to_up_sent_;
    int packets_to_down_sent_;
    int packets_to_left_sent_;
    int packets_to_right_sent_;

    std::vector<std::pair<int, int>> alltoall_send_matrix_;
    std::vector<std::pair<int, int>> alltoall_recv_matrix_;

    int alltoall_packet_sent_;
    int alltoall_packet_recved_;
};

}  // namespace AstraSim

#endif /* __MESHXY_HH__ */