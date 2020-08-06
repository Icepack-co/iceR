
#include "tsp-mcvfz472gty6.pb.h"
#include "tsptw-kcxbievqo879.pb.h"
#include "ivr7-kt461v8eoaif.pb.h"
#include "ivr8-yni1c9k2swof.pb.h"
#include "nvd-hap0j2y4zlm1.pb.h"
#include "matrix-vyv95n7wchpl.pb.h"
#include "ns3-tbfvuwtge2iq.pb.h"
#include "wrap.h"
#include <unordered_map>
using namespace std;

tspTabular tabulateTSPdata(string& tspSolveRequest, string& solRespString){
  TSP::SolveRequest r;

  if(r.ParseFromString(tspSolveRequest)){
    auto &m = r.model();

    tspTabular t;
    unordered_map<string, int> idlookup;
    for(int i = 0; i < m.points_size(); i++){
      auto &p =  m.points(i);
      idlookup.insert(pair<string,int>(p.id(), i));
    }

    TSP::SolutionResponse sol;
    if(sol.ParseFromString(solRespString)){
      t.nodes.resize(sol.tour_size());
      for(int i = 0; i < sol.tour_size(); i++){
        auto &tr = sol.tour(i);
        auto it = idlookup.find(tr);
        if(it != idlookup.end()){
          auto& p = m.points(it->second);
          t.nodes.sequence[i] = i+1;
          t.nodes.id[i] = p.id();
          t.nodes.x[i] = p.x();
          t.nodes.y[i]= p.y();
        }
      }
      t.edges.resize(sol.edges_size());
      for(int i =0; i < sol.edges_size(); i++){
        auto &e = sol.edges(i);
        geomvecs mygeom;
        mygeom.xs.resize(e.geometry_size());
        mygeom.ys.resize(e.geometry_size());
        for(int j =0; j <  e.geometry_size(); j++){
          mygeom.xs[j] = e.geometry(j).x();
          mygeom.ys[j] = e.geometry(j).y();
        }
        t.edges.geoms[i] = mygeom;
        t.edges.distance[i] = e.distance();
        t.edges.fromId[i] = e.from();
        t.edges.toId[i] = e.to();
      }
      return t;
    }
  }
  return tspTabular();
}

tspTabular tabulateTSPTWdata(string& tspSolveRequest, string& solRespString){
  TSPTW::SolveRequest r;

  if(r.ParseFromString(tspSolveRequest)){
    auto &m = r.model();

    tspTabular t;
    unordered_map<string, int> idlookup;
    for(int i = 0; i < m.points_size(); i++){
      auto &p =  m.points(i);
      idlookup.insert(pair<string,int>(p.id(), i));
    }

    TSPTW::SolutionResponse sol;
    if(sol.ParseFromString(solRespString)){
      t.nodes.resize(sol.tour_size());
      for(int i = 0; i < sol.tour_size(); i++){
        auto &tr = sol.tour(i);
        auto it = idlookup.find(tr);
        if(it != idlookup.end()){
          auto& p = m.points(it->second);
          t.nodes.sequence[i] = i+1;
          t.nodes.id[i] = p.id();
          t.nodes.x[i] = p.x();
          t.nodes.y[i]= p.y();
          t.nodes.arrivalTime[i] = sol.arrivaltimes(i);
          t.nodes.windowStart[i] = p.windowstart();
          t.nodes.windowEnd[i] = p.windowend();
        }
      }
      t.edges.resize(sol.edges_size());
      for(int i =0; i < sol.edges_size(); i++){
        auto &e = sol.edges(i);
        geomvecs mygeom;
        mygeom.xs.resize(e.geometry_size());
        mygeom.ys.resize(e.geometry_size());
        for(int j =0; j <  e.geometry_size(); j++){
          mygeom.xs[j] = e.geometry(j).x();
          mygeom.ys[j] = e.geometry(j).y();
        }
        t.edges.geoms[i] = mygeom;
        t.edges.distance[i] = e.distance();
        t.edges.time[i] = e.time();
        t.edges.fromId[i] = e.from();
        t.edges.toId[i] = e.to();
      }
      return t;
    }
  }
  return tspTabular();
}

ivr_tabular tabulateIVR7(string& ivrSolveRequest, string& solRespString){
  IVR7::SolutionResponse sr;
  if (sr.ParseFromString(solRespString)) {

    IVR7::SolveRequest r;
    if(r.ParseFromString(ivrSolveRequest)){
      unordered_map<string, int> idlookup;
      auto m = r.model();
      for(int i = 0; i < m.locations_size(); i++){
        auto &p =  m.locations(i);
        idlookup.insert(pair<string,int>(p.id(), i));
      }

      // first get all the dimension names
      unordered_set<string> dimensions;
      for (int i = 0; i < sr.routes_size(); i++) {
        auto& r = sr.routes(i);
        for (int j = 0; j < r.stops_size(); j++) {
          auto& s = r.stops(j);
          for (int k = 0; k < s.attributes_size(); k++) {
            dimensions.insert(s.attributes(k).dimid());
          }
        }
      }
      ivr_tabular t;
      // set up a row block with each attribute mask
      // run through all the route_rows. add em all up.
      for (int i = 0; i < sr.routes_size(); i++) {
        auto& r = sr.routes(i);
        for (int j = 0; j < r.stops_size(); j++) {
          auto& s = r.stops(j);
          stop_row rw;
          rw.stopid = s.id();
          rw.jobid = s.jobid();
          rw.location = s.locationid();
          rw.sequence = s.sequence();
          rw.taskid = s.taskid();
          rw.vehicleid = r.vehicleid();
          auto it = idlookup.find(rw.location);
          if(it != idlookup.end()){
            auto& p = m.locations(it->second);
            rw.x = p.geocode().longitude();
            rw.y= p.geocode().latitude();
          }

          for (auto dim : dimensions) {
            rw.dimVals.insert(pair<string, stop_block>(dim, stop_block()));
          }

          for (int k = 0; k < s.attributes_size(); k++) {
            auto& a = s.attributes(k);
            string d = a.dimid();
            auto& item = rw.dimVals[d];
            if (a.has_startvalue()) {
              item.start = a.startvalue();
            }
            if (a.has_endvalue()) {
              item.end = a.endvalue();
            }
            if (a.has_slackvalue()) {
              item.slackval = a.slackvalue();
            }
            if (a.has_slackcost()) {
              item.slackcost = a.slackcost();
            }
            if (a.has_tardyvalue()) {
              item.tardyval = a.tardyvalue();
            }
            if (a.has_tardycost()) {
              item.tardycost = a.tardycost();
            }
            if (a.has_cost()) {
              item.cost = a.cost();
            }
          }
          t.route_rows.push_back(rw);
        }

        for (int j = 0; j < r.interstops_size(); j++) {
          auto& is = r.interstops(j);
          edge_row er;
          er.fromStopId = is.fromstopid();
          er.toStopId = is.tostopid();
          er.geom.xs.resize(is.routesegments_size());
          er.geom.ys.resize(is.routesegments_size());
          for(int g = 0; g <  is.routesegments_size(); g++){
            auto& seg = is.routesegments(g);
            er.geom.xs[g] = seg.longitude();
            er.geom.ys[g] = seg.latitude();
          }

          for (auto dim : dimensions) {
            er.dimVals.insert(pair<string, edge_block>(dim, edge_block()));
          }
          for (int k = 0; k < is.attributes_size(); k++) {
            auto& a = is.attributes(k);
            string d = a.dimid();
            auto& item = er.dimVals[d];
            if (a.has_startvalue()) {
              item.start = a.startvalue();
            }
            if (a.has_endvalue()) {
              item.end = a.endvalue();
            }
            if (a.has_cost()) {
              item.cost = a.cost();
            }
          }
          t.edge_rows.push_back(er);
        }
      }
      return t;
    }
  }
  return ivr_tabular();
}


ivr_tabular tabulateIVR8(string& ivrSolveRequest, string& solRespString){
  IVR8::SolutionResponse sr;
  if (sr.ParseFromString(solRespString)) {

    IVR8::SolveRequest r;
    if(r.ParseFromString(ivrSolveRequest)){
      unordered_map<string, int> idlookup;
      auto m = r.model();
      for(int i = 0; i < m.locations_size(); i++){
        auto &p =  m.locations(i);
        idlookup.insert(pair<string,int>(p.id(), i));
      }

      // first get all the dimension names
      unordered_set<string> dimensions;
      for (int i = 0; i < sr.routes_size(); i++) {
        auto& r = sr.routes(i);
        for (int j = 0; j < r.stops_size(); j++) {
          auto& s = r.stops(j);
          for (int k = 0; k < s.attributes_size(); k++) {
            dimensions.insert(s.attributes(k).dimid());
          }
        }
      }
      ivr_tabular t;
      // set up a row block with each attribute mask
      // run through all the route_rows. add em all up.
      for (int i = 0; i < sr.routes_size(); i++) {
        auto& r = sr.routes(i);
        for (int j = 0; j < r.stops_size(); j++) {
          auto& s = r.stops(j);
          stop_row rw;
          rw.stopid = s.id();
          rw.jobid = s.jobid();
          rw.location = s.locationid();
          rw.sequence = s.sequence();
          rw.taskid = s.taskid();
          rw.vehicleid = r.vehicleid();
          rw.compartmentid = s.compartmentid();
          auto it = idlookup.find(rw.location);
          if(it != idlookup.end()){
            auto& p = m.locations(it->second);
            rw.x = p.geocode().longitude();
            rw.y= p.geocode().latitude();
          }

          for (auto dim : dimensions) {
            rw.dimVals.insert(pair<string, stop_block>(dim, stop_block()));
          }

          for (int k = 0; k < s.attributes_size(); k++) {
            auto& a = s.attributes(k);
            string d = a.dimid();
            auto& item = rw.dimVals[d];
            if (a.has_startvalue()) {
              item.start = a.startvalue();
            }
            if (a.has_endvalue()) {
              item.end = a.endvalue();
            }
            if (a.has_slackvalue()) {
              item.slackval = a.slackvalue();
            }
            if (a.has_slackcost()) {
              item.slackcost = a.slackcost();
            }
            if (a.has_tardyvalue()) {
              item.tardyval = a.tardyvalue();
            }
            if (a.has_tardycost()) {
              item.tardycost = a.tardycost();
            }
            if (a.has_cost()) {
              item.cost = a.cost();
            }
          }
          t.route_rows.push_back(rw);
        }

        for (int j = 0; j < r.interstops_size(); j++) {
          auto& is = r.interstops(j);
          edge_row er;
          er.fromStopId = is.fromstopid();
          er.toStopId = is.tostopid();
          er.geom.xs.resize(is.routesegments_size());
          er.geom.ys.resize(is.routesegments_size());
          for(int g = 0; g <  is.routesegments_size(); g++){
            auto& seg = is.routesegments(g);
            er.geom.xs[g] = seg.longitude();
            er.geom.ys[g] = seg.latitude();
          }

          for (auto dim : dimensions) {
            er.dimVals.insert(pair<string, edge_block>(dim, edge_block()));
          }
          for (int k = 0; k < is.attributes_size(); k++) {
            auto& a = is.attributes(k);
            string d = a.dimid();
            auto& item = er.dimVals[d];
            if (a.has_startvalue()) {
              item.start = a.startvalue();
            }
            if (a.has_endvalue()) {
              item.end = a.endvalue();
            }
            if (a.has_cost()) {
              item.cost = a.cost();
            }
          }
          t.edge_rows.push_back(er);
        }
      }
      return t;
    }
  }
  return ivr_tabular();
}

nvd_tabular tabulateNVD(string& nvdSolveRequest, string& nvdRespString){
  NVD::SolutionResponse sr;
  if (sr.ParseFromString(nvdRespString)) {
    nvd_tabular t;
    if(sr.has_instance()){
      auto soln = sr.instance();
      NVD::SolveRequest r;
      if(r.ParseFromString(nvdSolveRequest)){
        unordered_map<string, int> idlookup;
        auto m = r.model();
        for(int i = 0; i < m.visits_size(); i++){
          auto &p =  m.visits(i);
          idlookup.insert(pair<string,int>(p.id(), i));
        }
        for(int i =0 ; i < m.territories_size(); i++){
          auto &l = m.territories(i);
          idlookup.insert(pair<string,int>(l.id(), -i));
        }

        // first get all the dimension names
        unordered_set<string> dimensions;
        for (int i = 0; i < soln.routes_size(); i++) {
          auto& r = soln.routes(i);
          for (int j = 0; j < r.stops_size(); j++) {
            auto& s = r.stops(j);
            for (int k = 0; k < s.attributes_size(); k++) {
              dimensions.insert(s.attributes(k).dimid());
            }
          }
        }

        // set up a row block with each attribute mask
        // run through all the route_rows. add em all up.
        for (int i = 0; i < soln.routes_size(); i++) {
          auto& r = soln.routes(i);
          int dayIndex = r.day() + 1; // convert to 1-based for R.
          for (int j = 0; j < r.stops_size(); j++) {
            auto& s = r.stops(j);
            stop_row rw;
            rw.dayIndex = dayIndex;
            rw.stopid = s.id();
            rw.jobid = "";
            rw.location = s.locationid();
            rw.sequence = s.sequence();
            rw.taskid = s.visitid();
            rw.vehicleid = r.vehicleid();
            auto it = idlookup.find(rw.location);
            if(it != idlookup.end()){
              if(it->second > 0){
                auto& p = m.visits(it->second);
                rw.x = p.location().longitude();
                rw.y = p.location().latitude();
              }else{
                auto& l = m.territories(-it->second);
                rw.x = l.location().longitude();
                rw.y = l.location().latitude();
              }
            }

            for (auto dim : dimensions) {
              rw.dimVals.insert(pair<string, stop_block>(dim, stop_block()));
            }

            for (int k = 0; k < s.attributes_size(); k++) {
              auto& a = s.attributes(k);
              string d = a.dimid();
              auto& item = rw.dimVals[d];
              if (a.has_startvalue()) {
                item.start = a.startvalue();
              }
              if (a.has_endvalue()) {
                item.end = a.endvalue();
              }
              if (a.has_slackvalue()) {
                item.slackval = a.slackvalue();
              }
              if (a.has_slackcost()) {
                item.slackcost = a.slackcost();
              }
              if (a.has_tardyvalue()) {
                item.tardyval = a.tardyvalue();
              }
              if (a.has_tardycost()) {
                item.tardycost = a.tardycost();
              }
              if (a.has_cost()) {
                item.cost = a.cost();
              }
            }
            t.tab.route_rows.push_back(rw);
          }

          for (int j = 0; j < r.interstops_size(); j++) {
            auto& is = r.interstops(j);
            edge_row er;
            er.fromStopId = is.fromstopid();
            er.toStopId = is.tostopid();
            er.geom.xs.resize(is.routesegments_size());
            er.geom.ys.resize(is.routesegments_size());
            for(int g = 0; g <  is.routesegments_size(); g++){
              auto& seg = is.routesegments(g);
              er.geom.xs[g] = seg.longitude();
              er.geom.ys[g] = seg.latitude();
            }

            for (auto dim : dimensions) {
              er.dimVals.insert(pair<string, edge_block>(dim, edge_block()));
            }
            for (int k = 0; k < is.attributes_size(); k++) {
              auto& a = is.attributes(k);
              string d = a.dimid();
              auto& item = er.dimVals[d];
              if (a.has_startvalue()) {
                item.start = a.startvalue();
              }
              if (a.has_endvalue()) {
                item.end = a.endvalue();
              }
              if (a.has_cost()) {
                item.cost = a.cost();
              }
            }
            t.tab.edge_rows.push_back(er);
          }
        }
      }
    }
    if(sr.frontier_size() > 0){
      // then add all the frontier items.
      for(int i = 0; i < sr.frontier_size(); i++){
        auto f = sr.frontier(i);
        nvd_frontier_item fi;
        for(int j = 0; j < f.objectives_size(); j++){
          fi.objectiveValues.push_back(f.objectives(j));
          fi.objectiveNames.push_back(f.objectivenames(j));
        }
        fi.solutionIndex = i + 1;
        t.frontier.push_back(fi);
      }
    }
    return t;
  }
  return nvd_tabular();

}

matrix_tabular tabulateMatrix(string& matrixRequest, string& matrixResp){
  Matrix::MatrixRequest sr;
  if (sr.ParseFromString(matrixRequest)) {
    Matrix::MatrixResponse r;
    if(r.ParseFromString(matrixResp)){
      matrix_tabular t;
      for(int i =0 ; i < r.elements_size(); i++){
        auto &e = r.elements(i);
        t.fids.push_back(e.fromid());
        t.tids.push_back(e.toid());
        t.distances.push_back(e.distance());
        t.durations.push_back(e.duration());
      }
      return t;
    }
  }
  return matrix_tabular();
}

/// NS3
ns3_tabular tabulateNS3(string& ns3Request, string& ns3Resp){
  ns3_tabular n;
  NS3::SolveRequest sr;
  if(sr.ParseFromString(ns3Request)){
    auto& m = sr.model();
    for(int i = 0; i < m.nodes_size(); i++){
      auto& node = m.nodes(i);
      n.nodes.Id.push_back(node.id());
      n.nodes.xs.push_back(node.geocode().longitude());
      n.nodes.ys.push_back(node.geocode().latitude());
    }

    NS3::SolutionResponse r;
    if(r.ParseFromString(ns3Resp)){
      auto &na = n.assignmenttab;
      for(int i = 0; i < r.assignments_size(); i++){
        auto &a = r.assignments(i);
        na.source.push_back(a.source());
        na.destination.push_back(a.destination());
        na.productId.push_back(a.productid());
        na.amount.push_back(a.amount());
        na.cost.push_back(a.cost());
        na.laneRateId.push_back(a.lanerateid());
        na.costModelId.push_back(a.costmodelid());
        na.distance.push_back(a.distance());
        na.duration.push_back(a.duration());
      }
      auto &nfpt = n.nodeflowproducttab;
      for(int i = 0; i < r.nodeproductflows_size(); i++){
        auto &nf = r.nodeproductflows(i);
        nfpt.nodeId.push_back(nf.nodeid());
        nfpt.productId.push_back(nf.productid());
        nfpt.inFlow.push_back(nf.inflow());
        nfpt.outFlow.push_back(nf.outflow());
        nfpt.flowCost.push_back(nf.flowcost());
        nfpt.fixedCost.push_back(nf.fixedcost());
        nfpt.productionAmount.push_back(nf.productionamount());
        nfpt.productionPenalty.push_back(nf.productionpenalty());
        nfpt.productionCost.push_back(nf.productioncost());
        nfpt.consumptionAmount.push_back(nf.consumptionamount());
        nfpt.consumptionPenalty.push_back(nf.consumptionpenalty());
        nfpt.consumptionCost.push_back(nf.consumptioncost());
      }
      auto &nft = n.nodeflowtab;
      for(int i = 0; i < r.nodeflows_size(); i++){
        auto &nf = r.nodeflows(i);
        nft.nodeId.push_back(nf.nodeid());
        nft.inFlow.push_back(nf.inflow());
        nft.outFlow.push_back(nf.outflow());
        nft.flowCost.push_back(nf.flowcost());
        nft.fixedCost.push_back(nf.fixedcost());
        nft.productFlowCost.push_back(nf.productflowcost());
        nft.productFixedCost.push_back(nf.productfixedcost());
        nft.productionAmount.push_back(nf.productionamount());
        nft.productionPenalty.push_back(nf.productionpenalty());
        nft.productionCost.push_back(nf.productioncost());
        nft.consumptionAmount.push_back(nf.consumptionamount());
        nft.consumptionPenalty.push_back(nf.consumptionpenalty());
        nft.consumptionCost.push_back(nf.consumptioncost());
      }
      for(int i =0 ; i < r.routes_size(); i++){
        auto rt = r.routes(i);
        n.geomtab.fromId.push_back(rt.fromid());
        n.geomtab.toId.push_back(rt.toid());
        geomvecs gv;
        for(int j = 0; j < rt.geometrysequence_size(); j++){
          int gind = rt.geometrysequence(j);
          auto& gs = r.geometrysequence(gind);
          for(int p = 0; p < gs.x_size(); p++){
            gv.xs.push_back(gs.x(p));
            gv.ys.push_back(gs.y(p));
          }
        }
        n.geomtab.geom.push_back(gv);
      }
    }
  }
  return n;
}

