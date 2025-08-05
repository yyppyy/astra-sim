/******************************************************************************
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
*******************************************************************************/

#ifndef __RECV_PACKET_EVENT_HANDLER_DATA_HH__
#define __RECV_PACKET_EVENT_HANDLER_DATA_HH__

#include "astra-sim/system/BaseStream.hh"
#include "astra-sim/system/BasicEventHandlerData.hh"
#include "astra-sim/system/astraccl/custom_collectives/CustomAlgorithm.hh"

namespace AstraSim {

class WorkloadLayerHandlerData;

class RecvPacketEventHandlerData : public BasicEventHandlerData {
  public:
    RecvPacketEventHandlerData();
    RecvPacketEventHandlerData(BaseStream* owner,
                               int sys_id,
                               EventType event,
                               int vnet,
                               int stream_id);
    RecvPacketEventHandlerData(BaseStream* owner,
                               int sys_id,
                               EventType event,
                               int vnet,
                               int stream_id,
                               int tag);

    Workload* workload;
    WorkloadLayerHandlerData* wlhd;
    BaseStream* owner;
    CustomAlgorithm* custom_algorithm;
    int vnet;
    int stream_id;
    bool message_end;
    Tick ready_time;
    int tag;
};

}  // namespace AstraSim

#endif /* __RECV_PACKET_EVENT_HANDLER_DATA_HH__ */
