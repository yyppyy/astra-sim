/******************************************************************************
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
*******************************************************************************/

#ifndef __MESH_TOPOLOGY_HH__
#define __MESH_TOPOLOGY_HH__

#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "astra-sim/system/astraccl/native_collectives/logical_topology/BasicLogicalTopology.hh"

namespace AstraSim {

class MeshTopology : public BasicLogicalTopology {
  public:
    // enum class Direction { Clockwise, Anticlockwise };
    // enum class Dimension { Local, Vertical, Horizontal, NA };
    int get_num_of_nodes_in_dimension(int dimension) override;
    // RingTopology(Dimension dimension,
    //              int id,
    //              int total_nodes_in_ring,
    //              int index_in_ring,
    //              int offset);
    MeshTopology(int dimension, int id, int npu_count);
    // RingTopology(Dimension dimension, int id, std::vector<int> NPUs);
    // virtual int get_receiver(int node_id, Direction direction);
    // virtual int get_sender(int node_id, Direction direction);
    // int get_nodes_in_ring();
    // bool is_enabled();
    // Dimension get_dimension();
    // int get_index_in_ring();
    int get_x();
    int get_y();

  private:

    int x_;
    int y_;

    // std::unordered_map<int, int> id_to_index;
    // std::unordered_map<int, int> index_to_id;

    // std::string name;
    // int id;
    // int offset;
    // int total_nodes_in_ring;
    // int index_in_ring;
    // Dimension dimension;

    // virtual int get_receiver_homogeneous(int node_id,
    //                                      Direction direction,
    //                                      int offset);
};

}  // namespace AstraSim

#endif /* __MESH_TOPOLOGY_HH__ */
