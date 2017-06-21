/**
 * @file   GeometryIteratorLoopTestAlg.cxx
 * @brief  Tests the correct iteration of the geo::Geometry iterators
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   August 25, 2014
 */

// our header
#include "GeometryIteratorLoopTestAlg.h"

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"

// Framework includes
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace geo {

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

  //......................................................................
  unsigned int GeometryIteratorLoopTestAlg::Run() {
    
    /*
     * This monstrous test looks through the elements of the geometry.
     * It is structured as follows:
     * 
     * - loop on cryostats by index
     *   * check global cryostat iterators (ID and element)
     *   - loop on TPCs by index
     *     * check global TPC iterators (ID and element)
     *     * check local TPC iterators (cryostat)
     *     - loop on planes by index
     *       * check global plane iterators (ID and element)
     *       * check local plane iterators (cryostat)
     *       - loop on wires by index
     *         * check global wire iterators (ID and element)
     *         * check local wire iterators (cryostat)
     *         * increase wire iterators (including cryostat locals)
     *       * increase plane iterators (including cryostat locals)
     *     * increase TPC iterators (including cryostat locals)
     *   * check cryostat-local iterators (TPC, plane, wire IDs and elements)
     *   - loop on TPC sets by index
     *     * check global TPC set iterators (ID)
     *     - loop on readout planes by index
     *       * check global readout plane iterators (ID)
     *       - loop on channels by channel ID (currently disabled)
     *       * increase readout plane iterators
     *     * increase TPC set iterators
     * - loops by range-for
     *   * by cryostat ID
     *   * by cryostat
     *   * by TPC ID
     *   * by TPC
     *   * by plane ID
     *   * by plane
     *   * by wire ID
     *   * by wire
     *   * by TPC set ID
     *   * by TPC set
     *   * by readout plane ID
     *   * by readout plane
     * 
     * In words: the test is structured in two almost-independent parts.
     * In the first, nested loops are driven by element indices.
     * In the second, range-for loops are implemented. The number of loops in
     * here is compared to the number of iterations recorded in the first
     * section (hence the dependence of the second part from the first one).
     * For the elements with both ID and geometry class (cryostat, TPC, plane
     * and wire) both types of iterators, on ID and on element, are checked,
     * while the ones lacking an element class (TPC set, readout plane) only
     * the ID iterators are tested. The same holds for range-for loops too.
     * 
     * In each index loop, a loop of the contained element is nested. Also,
     * the iterators concerning the indexed element are checked. Finally,
     * range-for loops are rolled for iterators local to the element.
     * For example, each iteration of the TPC loop includes a wire plane loop,
     * a check on TPC iterators (both global and local to the cryostat), and
     * a range-for loop on planes and wires on the TPC.
     * 
     * Cryostat loop contains tests for both TPCs and TPC sets.
     * 
     */
    
    
    const unsigned int nCryo = geom->Ncryostats(); 
    LOG_VERBATIM("GeometryIteratorLoopTest")
      << "We have " << nCryo << " cryostats";
    
    unsigned int nErrors = 0;
    unsigned int nCryostats = 0, nTPCs = 0, nPlanes = 0, nWires = 0,
      cumTPCsets = 0, cumROPs = 0;
    geo::GeometryCore::cryostat_id_iterator iCryostatID
      = geom->begin_cryostat_id();
    geo::GeometryCore::TPC_id_iterator iTPCID = geom->begin_TPC_id();
    geo::GeometryCore::plane_id_iterator iPlaneID = geom->begin_plane_id();
    geo::GeometryCore::wire_id_iterator iWireID = geom->begin_wire_id();
    
    auto iTPCsetID = geom->begin_TPCset_id();
    auto iROPID = geom->begin_ROP_id();
    // no channel iterator implemented
    
    geo::GeometryCore::cryostat_iterator iCryostat = geom->begin_cryostat();
    geo::GeometryCore::TPC_iterator iTPC = geom->begin_TPC();
    geo::GeometryCore::plane_iterator iPlane = geom->begin_plane();
    geo::GeometryCore::wire_iterator iWire = geom->begin_wire();
    
    geo::CryostatID runningCID{ 0 };
    geo::TPCID runningTID{ runningCID, 0 };
    geo::PlaneID runningPID{ runningTID, 0 };
    geo::WireID runningWID{ runningPID, 0 };
    readout::TPCsetID runningSID{ runningCID, 0 };
    readout::ROPID runningRID{ runningSID, 0 };
    
    for(unsigned int c = 0; c < nCryo; ++c) {
      const geo::CryostatID expCID(c);
      const CryostatGeo& cryo(geom->Cryostat(c));
      const unsigned int nTPC = cryo.NTPC();
      const unsigned int nTPCsets = geom->NTPCsets(expCID);
      
      LOG_TRACE("GeometryIteratorLoopTest") << "  C=" << c
        << " (" << nTPC << " TPCs, " << nTPCsets << " TPC sets)";
      
      if (runningCID != expCID) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "Cryostat ID incremented to " << runningCID << ", expected: "
          << expCID;
        ++nErrors;
      }
      
      if (!iCryostatID) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "Cryostat ID iterator thinks it's all over at C=" << c;
        ++nErrors;
      }
      else if (iCryostatID->Cryostat != c) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "Cryostat ID iterator thinks it's at C=" << (*iCryostatID)
          << " instead of " << c;
        ++nErrors;
      }
      else if (iCryostatID.get() != &cryo) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "Cryostat ID iterator retrieves CryostatGeo["
          << ((void*) iCryostatID.get())
          << "] instead of [" << ((void*) &cryo) << "]";
        ++nErrors;
      }
      
      if (&*iCryostat != &cryo) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "Cryostat iterator retrieves CryostatGeo["
          << ((void*) iCryostat.get())
          << "] (" << iCryostat.ID() << ") instead of ["
          << ((void*) &cryo) << "] (C=" << c << ")";
        ++nErrors;
      }
      
      geo::TPC_id_iterator iTPCIDinCryo = geom->begin_TPC_id(expCID);
      geo::TPC_iterator iTPCinCryo = geom->begin_TPC(expCID);
      geo::plane_id_iterator iPlaneIDinCryo = geom->begin_plane_id(expCID);
      geo::plane_iterator iPlaneInCryo = geom->begin_plane(expCID);
      geo::wire_id_iterator iWireIDinCryo = geom->begin_wire_id(expCID);
      geo::wire_iterator iWireInCryo = geom->begin_wire(expCID);
      
      unsigned int nPlanesInCryo = 0;
      unsigned int nWiresInCryo = 0;
      
      for(unsigned int t = 0; t < nTPC; ++t){
        const TPCGeo& TPC(cryo.TPC(t));
        const geo::TPCID expTID(expCID, t);
        const unsigned int NPlanes = TPC.Nplanes();
        
        LOG_TRACE("GeometryIteratorLoopTest") << "    " << expTID
          << " (" << NPlanes << " planes)";
        
        if (runningTID != expTID) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "TPC ID incremented to " << runningTID << ", expected: "
            << expTID;
          ++nErrors;
        }
        
        if (!iTPCID) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "TPC ID iterator thinks it's all over at " << expTID;
          ++nErrors;
        }
        else if (iTPCID->Cryostat != c) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "TPC ID iterator thinks it's at C=" << iTPCID->Cryostat
            << " instead of " << c;
          ++nErrors;
        }
        else if (iTPCID->TPC != t) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "TPC ID iterator thinks it's at T=" << iTPCID->TPC
            << " instead of " << t;
          ++nErrors;
        }
        else if (iTPCID.get() != &TPC) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "TPC ID iterator retrieves TPCGeo[" << ((void*) iTPCID.get())
            << "] instead of [" << ((void*) &TPC) << "]";
          ++nErrors;
        }
        
        if (&*iTPC != &TPC) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "TPC iterator retrieves TPCGeo[" << ((void*) iTPC.get())
            << "] (" << iTPC.ID() << ") instead of [" << ((void*) &TPC)
            << "] (" << expTID << ")";
          ++nErrors;
        }
        
        if (*iTPCIDinCryo != expTID) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "TPC ID local iterator in " << expCID
            << " points to " << *iTPCIDinCryo << " instead of " << expTID;
          ++nErrors;
        }
        if (iTPCinCryo->ID() != expTID) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "TPC local iterator in " << expCID
            << " points to " << iTPCinCryo->ID() << " instead of " << expTID;
          ++nErrors;
        }
        
        for(unsigned int p = 0; p < NPlanes; ++p) {
          const PlaneGeo& Plane(TPC.Plane(p));
          geo::PlaneID const expPID(expTID, p);
          const unsigned int NWires = Plane.Nwires();
          
          LOG_TRACE("GeometryIteratorLoopTest") << "    " << expTID
            << " (" << NWires << " wires)";
          
          if (runningPID != expPID) {
            LOG_ERROR("GeometryIteratorLoopTest")
              << "Plane ID incremented to " << runningPID << ", expected: "
              << expPID;
            ++nErrors;
          }
          if (!iPlaneID) {
            LOG_ERROR("GeometryIteratorLoopTest")
              << "plane ID iterator thinks it's all over at " << expPID;
            ++nErrors;
          }
          else if (iPlaneID->Cryostat != c) {
            LOG_ERROR("GeometryIteratorLoopTest")
              << "plane ID iterator thinks it's at C=" << iPlaneID->Cryostat
              << " instead of " << c;
            ++nErrors;
          }
          else if (iPlaneID->TPC != t) {
            LOG_ERROR("GeometryIteratorLoopTest")
              << "plane ID iterator thinks it's at T=" << iPlaneID->TPC
              << " instead of " << t;
            ++nErrors;
          }
          else if (iPlaneID->Plane != p) {
            LOG_ERROR("GeometryIteratorLoopTest")
              << "plane ID iterator thinks it's at P=" << iPlaneID->Plane
              << " instead of " << p;
            ++nErrors;
          }
          else if (iPlaneID.get() != &Plane) {
            LOG_ERROR("GeometryIteratorLoopTest")
              << "plane ID iterator retrieves PlaneGeo["
              << ((void*) iPlaneID.get()) << "] instead of ["
               << ((void*) &Plane) << "]";
            ++nErrors;
          }
          
          if (&*iPlane != &Plane) {
            LOG_ERROR("GeometryIteratorLoopTest")
              << "plane iterator retrieves PlaneGeo[" << ((void*) iPlane.get())
              << "] instead of [" << ((void*) &Plane) << "] (" << expPID << ")";
            ++nErrors;
          }
          
          if (*iPlaneIDinCryo != expPID) {
            LOG_ERROR("GeometryIteratorLoopTest")
              << "Plane ID local iterator in " << expCID << " points to "
              << *iPlaneIDinCryo << " instead of " << expPID;
            ++nErrors;
          }
          if (iPlaneInCryo->ID() != expPID) {
            LOG_ERROR("GeometryIteratorLoopTest")
              << "Plane local iterator in " << expCID << " points to "
              << iPlaneInCryo->ID() << " instead of " << expPID;
            ++nErrors;
          }
          
          
          for(unsigned int w = 0; w < NWires; ++w) {
            const WireGeo& Wire(Plane.Wire(w));
            geo::WireID const expWID(expPID, w);
            
            if (runningWID != expWID) {
              LOG_ERROR("GeometryIteratorLoopTest")
                << "Wire ID incremented to " << runningWID << ", expected: "
                << expWID;
              ++nErrors;
            }
            
            LOG_TRACE("GeometryIteratorLoopTest") << "    " << expWID;
            if (!iWireID) {
              LOG_ERROR("GeometryIteratorLoopTest")
                << "wire ID iterator thinks it's all over at " << expWID;
              ++nErrors;
            }
            else if (iWireID->Cryostat != c) {
              LOG_ERROR("GeometryIteratorLoopTest")
                << "wire ID iterator thinks it's at C=" << iWireID->Cryostat
                << " instead of " << c;
              ++nErrors;
            }
            else if (iWireID->TPC != t) {
              LOG_ERROR("GeometryIteratorLoopTest")
                << "wire ID iterator thinks it's at T=" << iWireID->TPC
                << " instead of " << t;
              ++nErrors;
            }
            else if (iWireID->Plane != p) {
              LOG_ERROR("GeometryIteratorLoopTest")
                << "wire ID iterator thinks it's at P=" << iWireID->Plane
                << " instead of " << p;
              ++nErrors;
            }
            else if (iWireID->Wire != w) {
              LOG_ERROR("GeometryIteratorLoopTest")
                << "wire ID iterator thinks it's at W=" << iWireID->Wire
                << " instead of " << w;
              ++nErrors;
            }
            else if (iWireID.get() != &Wire) {
              LOG_ERROR("GeometryIteratorLoopTest")
                << "wire ID iterator retrieves WireGeo["
                << ((void*) iWireID.get())
                << "] instead of [" << ((void*) &Wire) << "]";
              ++nErrors;
            }
            
            if (&*iWire != &Wire) {
              LOG_ERROR("GeometryIteratorLoopTest")
                << "wire iterator retrieves WireGeo[" << ((void*) iWire.get())
                << "] instead of [" << ((void*) &Wire) << "] (" << expWID
                << ")";
              ++nErrors;
            }
          
            if (*iWireIDinCryo != expWID) {
              LOG_ERROR("GeometryIteratorLoopTest")
                << "Wire ID local iterator in " << expCID << " points to "
                << *iWireIDinCryo << " instead of " << expWID;
              ++nErrors;
            }
            if (iWireInCryo.ID() != expWID) {
              LOG_ERROR("GeometryIteratorLoopTest")
                << "Wire local iterator in " << expCID << " points to "
                << iWireInCryo.ID() << " instead of " << expWID;
              ++nErrors;
            }
          
            ++iWireID;
            ++iWire;
            ++iWireIDinCryo;
            ++iWireInCryo;
            ++nWires;
            ++nWiresInCryo;
            geom->IncrementID(runningWID);
          } // end loop over wires
          ++iPlaneID;
          ++iPlane;
          ++iPlaneIDinCryo;
          ++iPlaneInCryo;
          ++nPlanes;
          ++nPlanesInCryo;
          geom->IncrementID(runningPID);
        } // end loop over planes
        
        ++iTPCID;
        ++iTPC;
        ++iTPCIDinCryo;
        ++iTPCinCryo;
        ++nTPCs;
        geom->IncrementID(runningTID);
      } // end loop over tpcs
      
      if (iTPCIDinCryo != geom->end_TPC_id(expCID)) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "TPC ID local iterator in " << expCID
          << " should be at end, and instead points to " << *iTPCIDinCryo;
        ++nErrors;
      }
      if (iTPCinCryo != geom->end_TPC(expCID)) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "TPC local iterator in " << expCID
          << " should be at end, and instead points to " << iTPCinCryo->ID();
        ++nErrors;
      }
      
      if (iPlaneIDinCryo != geom->end_plane_id(expCID)) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "Plane ID local iterator in " << expCID
          << " should be at end, and instead points to " << *iPlaneIDinCryo;
        ++nErrors;
      }
      if (iPlaneInCryo != geom->end_plane(expCID)) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "Plane local iterator in " << expCID
          << " should be at end, and instead points to " << iPlaneInCryo->ID();
        ++nErrors;
      }
      
      if (iWireIDinCryo != geom->end_wire_id(expCID)) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "Wire ID local iterator in " << expCID
          << " should be at end, and instead points to " << *iWireIDinCryo;
        ++nErrors;
      }
      if (iWireInCryo != geom->end_wire(expCID)) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "Wire local iterator in " << expCID
          << " should be at end, and instead points to " << iWireInCryo.ID();
        ++nErrors;
      }
      
      
      for(unsigned int s = 0; s < nTPCsets; ++s) {
        readout::TPCsetID const expSID(expCID, s);
        const unsigned int NROPs = geom->NROPs(expSID);
        
        LOG_TRACE("GeometryIteratorLoopTest")
          << "    " << expSID << " (" << NROPs << " planes)";
        
        if (runningSID != expSID) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "TPC set ID incremented to " << runningSID << ", expected: "
            << expSID;
          ++nErrors;
        }
        
        if (!iTPCsetID) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "TPCset ID iterator thinks it's all over at " << expSID;
          ++nErrors;
        }
        else if (iTPCsetID->Cryostat != c) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "TPC set ID iterator thinks it's at C=" << iTPCsetID->Cryostat
            << " instead of " << c;
          ++nErrors;
        }
        else if (iTPCsetID->TPCset != s) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "TPC set ID iterator thinks it's at S=" << iTPCsetID->TPCset
            << " instead of " << s;
          ++nErrors;
        }
        
        for(unsigned int r = 0; r < NROPs; ++r) {
          readout::ROPID const expRID(expSID, r);
          const unsigned int NChannels = geom->Nchannels(expRID);
          
          LOG_TRACE("GeometryIteratorLoopTest")
            << "    " << expRID << " (" << NChannels << " channels)";
          
          if (runningRID != expRID) {
            LOG_ERROR("GeometryIteratorLoopTest")
              << "Readout plane ID incremented to " << runningRID
              << ", expected: " << expRID;
            ++nErrors;
          }
          
          if (!iROPID) {
            LOG_ERROR("GeometryIteratorLoopTest")
              << "readout plane ID iterator thinks it's all over at " << expRID;
            ++nErrors;
          }
          else if (iROPID->Cryostat != c) {
            LOG_ERROR("GeometryIteratorLoopTest")
              << "readout plane ID iterator thinks it's at C="
              << iROPID->Cryostat << " instead of " << c;
            ++nErrors;
          }
          else if (iROPID->TPCset != s) {
            LOG_ERROR("GeometryIteratorLoopTest")
              << "readout plane ID iterator thinks it's at S=" << iROPID->TPCset
              << " instead of " << s;
            ++nErrors;
          }
          else if (iROPID->ROP != r) {
            LOG_ERROR("GeometryIteratorLoopTest")
              << "readout plane ID iterator thinks it's at R=" << iROPID->ROP
              << " instead of " << r;
            ++nErrors;
          }
          
          /*
          for(unsigned int ch = 0; ch < NChannels; ++ch) {
            
            LOG_TRACE("GeometryIteratorLoopTest") << "    channel=" << ch;
          
            ++iChannelID;
            ++nChannels;
          } // end loop over channels
          */
          
          ++iROPID;
          ++cumROPs;
          geom->IncrementID(runningRID);
        } // end loop over readout planes
        ++iTPCsetID;
        ++cumTPCsets;
        geom->IncrementID(runningSID);
      } // end loop over TPC sets
      
      //
      // test if we can loop all TPCs in this cryostat via iterator box
      //
      unsigned int nTPCsInCryo = cryo.NTPC();
      LOG_DEBUG("GeometryIteratorsDump")
        << "Looping though " << nTPCsInCryo << " TPCs in " << expCID;
      unsigned int nLoopedTPCIDs = 0;
      for (geo::TPCID const& tID: geom->IterateTPCIDs(expCID)) {
        LOG_TRACE("GeometryIteratorsDump") << tID;
        if (nLoopedTPCIDs >= nTPCsInCryo) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "After all " << nLoopedTPCIDs << " TPC IDs in " << expCID
            << ", iterator has not reached the end ("
            << *(geom->end_TPC_id(expCID)) << ") but it's still at " << tID;
          ++nErrors;
          break;
        }
        ++nLoopedTPCIDs;
      }
      if (nLoopedTPCIDs < nTPCsInCryo) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "Looped only " << nLoopedTPCIDs
          << " TPC IDs in " << expCID << ", while we expected "
          << nTPCsInCryo << " iterations!";
        ++nErrors;
      } // if
      unsigned int nLoopedTPCs = 0;
      for (geo::TPCGeo const& TPC: geom->IterateTPCs(expCID)) {
        if (nLoopedTPCs >= nTPCsInCryo) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "After all " << nLoopedTPCs
            << " TPCs in " << expCID << ", iterator has not reached the end";
          ++nErrors;
          break;
        }
        ++nLoopedTPCs;
      }
      if (nLoopedTPCs < nTPCsInCryo) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "Looped only " << nLoopedTPCs << " TPCs in " << expCID
          << ", while we expected " << nTPCsInCryo << " iterations!";
        ++nErrors;
      } // if
      
      //
      // test if we can loop all planes in this cryostat via iterator box
      //
      LOG_DEBUG("GeometryIteratorsDump")
        << "Looping though " << nPlanesInCryo << " planes in " << expCID;
      unsigned int nLoopedPlaneIDs = 0;
      for (geo::PlaneID const& pID: geom->IteratePlaneIDs(expCID)) {
        LOG_TRACE("GeometryIteratorsDump") << pID;
        if (nLoopedPlaneIDs >= nPlanesInCryo) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "After all " << nLoopedPlaneIDs << " plane IDs in " << expCID
            << ", iterator has not reached the end ("
            << *(geom->end_plane_id(expCID)) << ") but it's still at " << pID;
          ++nErrors;
          break;
        }
        ++nLoopedPlaneIDs;
      }
      if (nLoopedPlaneIDs < nPlanesInCryo) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "Looped only " << nLoopedPlaneIDs
          << " plane IDs in " << expCID << ", while we expected "
          << nPlanesInCryo << " iterations!";
        ++nErrors;
      } // if
      unsigned int nLoopedPlanes = 0;
      for (geo::PlaneGeo const& plane: geom->IteratePlanes(expCID)) {
        if (nLoopedPlanes >= nPlanesInCryo) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "After all " << nLoopedPlanes
            << " planes in " << expCID << ", iterator has not reached the end";
          ++nErrors;
          break;
        }
        ++nLoopedPlanes;
      }
      if (nLoopedPlanes < nPlanesInCryo) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "Looped only " << nLoopedPlanes << " planes in " << expCID
          << ", while we expected " << nPlanesInCryo << " iterations!";
        ++nErrors;
      } // if
      
      //
      // test if we can loop all wires in this cryostat via iterator box
      //
      LOG_DEBUG("GeometryIteratorsDump")
        << "Looping though " << nWiresInCryo << " wires in " << expCID;
      unsigned int nLoopedWireIDs = 0;
      for (geo::WireID const& wID: geom->IterateWireIDs(expCID)) {
        LOG_TRACE("GeometryIteratorsDump") << wID;
        if (nLoopedWireIDs >= nWiresInCryo) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "After all " << nLoopedWireIDs << " wire IDs in " << expCID
            << ", iterator has not reached the end ("
            << *(geom->end_wire_id(expCID)) << ") but it's still at " << wID;
          ++nErrors;
          break;
        }
        ++nLoopedWireIDs;
      }
      if (nLoopedWireIDs < nWiresInCryo) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "Looped only " << nLoopedWireIDs
          << " wire IDs in " << expCID << ", while we expected "
          << nWiresInCryo << " iterations!";
        ++nErrors;
      } // if
      unsigned int nLoopedWires = 0;
      for (geo::WireGeo const& wire: geom->IterateWires(expCID)) {
        if (nLoopedWires >= nWiresInCryo) {
          LOG_ERROR("GeometryIteratorLoopTest")
            << "After all " << nLoopedWires
            << " wires in " << expCID << ", iterator has not reached the end";
          ++nErrors;
          break;
        }
        ++nLoopedWires;
      }
      if (nLoopedWires < nWiresInCryo) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "Looped only " << nLoopedWires << " wires in " << expCID
          << ", while we expected " << nWiresInCryo << " iterations!";
        ++nErrors;
      } // if
      
      
      ++iCryostatID;
      ++iCryostat;
      ++nCryostats;
      geom->IncrementID(runningCID);
    } // end loop over cryostats
    
    if (runningCID) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "cryostat ID still valid (" << runningCID
        << ") after incrementing from the last one.";
      ++nErrors;
    }
    
    if (iCryostatID) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "Cryostat ID iterator thinks it's still at " << *iCryostatID
        << ", but we are already finished";
      ++nErrors;
    }
    try {
      geo::CryostatGeo const& Cryo = *iCryostat;
      LOG_ERROR("GeometryIteratorLoopTest")
        << "cryostat iterator thinks it's still at " << iCryostat.ID()
        << ", but we are already finished";
      ++nErrors;
    }
    catch (cet::exception const&) {
      LOG_DEBUG("GeometryIteratorLoopTest") << "exception caught"
        " while dereferencing an iterator to a past-the-end cryostat.\n";
    }
    
    // test if we can loop all cryostats with the iterators (via iterator box)
    LOG_DEBUG("GeometryIteratorsDump")
      << "Looping though " << nCryostats << " cryostats";
    unsigned int nLoopedCryostatIDs = 0;
    for (geo::CryostatID const& cID: geom->IterateCryostatIDs()) {
      LOG_TRACE("GeometryIteratorsDump") << cID;
      if (nLoopedCryostatIDs >= nCryostats) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "After all " << nLoopedCryostatIDs
          << " cryostat IDs, iterator has not reached the end ("
          << *(geom->end_cryostat_id()) << ") but it's still at " << cID;
        ++nErrors;
        break;
      }
      ++nLoopedCryostatIDs;
    }
    if (nLoopedCryostatIDs < nCryostats) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "Looped only " << nLoopedCryostatIDs
        << " cryostat IDs, while we expected " << nCryostats << " iterations!";
      ++nErrors;
    } // if
    unsigned int nLoopedCryostats = 0;
    for (geo::CryostatGeo const& Cryo: geom->IterateCryostats()) {
      if (nLoopedCryostats >= nCryostats) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "After all " << nLoopedCryostats
          << " cryostats, iterator has not reached the end";
        ++nErrors;
        break;
      }
      ++nLoopedCryostats;
    }
    if (nLoopedCryostats < nCryostats) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "Looped only " << nLoopedCryostats
        << " cryostats, while we expected " << nCryostats << " iterations!";
      ++nErrors;
    } // if
    
    if (runningTID) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "TPC ID still valid (" << runningTID
        << ") after incrementing from the last one.";
      ++nErrors;
    }
    
    if (iTPCID) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "TPC iterator thinks it's still at " << *iTPCID
        << ", but we are already finished";
      ++nErrors;
    }
    try {
      geo::TPCGeo const& TPC = *iTPC;
      LOG_ERROR("GeometryIteratorLoopTest")
        << "TPC iterator thinks it's still at " << iTPC.ID()
        << ", but we are already finished";
      ++nErrors;
    }
    catch (cet::exception const&) {
      LOG_DEBUG("GeometryIteratorLoopTest") << "exception caught"
        " while dereferencing an iterator to a past-the-end TPC.\n";
    }
    
    // test if we can loop all TPCs with the iterators (via iterator box)
    LOG_DEBUG("GeometryIteratorsDump")
      << "Looping though " << nTPCs << " TPCs";
    unsigned int nLoopedTPCIDs = 0;
    for (geo::TPCID const& tID: geom->IterateTPCIDs()) {
      LOG_TRACE("GeometryIteratorsDump") << tID;
      if (nLoopedTPCIDs >= nTPCs) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "After all " << nLoopedTPCIDs
          << " TPC IDs, iterator has not reached the end ("
          << *(geom->end_TPC_id()) << ") but it's still at " << tID;
        ++nErrors;
        break;
      }
      ++nLoopedTPCIDs;
    }
    if (nLoopedTPCIDs < nTPCs) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "Looped only " << nLoopedTPCIDs
        << " TPC IDs, while we expected " << nTPCs << " iterations!";
      ++nErrors;
    } // if
    unsigned int nLoopedTPCs = 0;
    for (geo::TPCGeo const& TPC: geom->IterateTPCs()) {
      if (nLoopedTPCs >= nTPCs) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "After all " << nLoopedTPCs
          << " TPCs, iterator has not reached the end";
        ++nErrors;
        break;
      }
      ++nLoopedTPCs;
    }
    if (nLoopedTPCs < nTPCs) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "Looped only " << nLoopedTPCs
        << " TPCs, while we expected " << nTPCs << " iterations!";
      ++nErrors;
    } // if
    
    
    if (runningPID) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "Plane ID still valid (" << runningPID
        << ") after incrementing from the last one.";
      ++nErrors;
    }
    
    if (iPlaneID) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "plane iterator thinks it's still at " << *iPlaneID
        << ", but we are already finished";
      ++nErrors;
    }
    try {
      geo::PlaneGeo const& Plane = *iPlane;
      LOG_ERROR("GeometryIteratorLoopTest")
        << "plane iterator thinks it's still at " << iPlane.ID()
        << ", but we are already finished";
      ++nErrors;
    }
    catch (cet::exception const&) {
      LOG_DEBUG("GeometryIteratorLoopTest") << "exception caught"
        " while dereferencing an iterator to a past-the-end plane.\n";
    }
    
    // test if we can loop all planes with the iterators (via iterator box)
    LOG_DEBUG("GeometryIteratorsDump")
      << "Looping though " << nPlanes << " planes";
    unsigned int nLoopedPlaneIDs = 0;
    for (geo::PlaneID const& pID: geom->IteratePlaneIDs()) {
      LOG_TRACE("GeometryIteratorsDump") << pID;
      if (nLoopedPlaneIDs >= nPlanes) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "After all " << nLoopedPlaneIDs
          << " planes, ID iterator has not reached the end ("
          << *(geom->end_plane_id()) << ") but it's still at " << pID;
        ++nErrors;
        break;
      }
      ++nLoopedPlaneIDs;
    }
    if (nLoopedPlaneIDs < nPlanes) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "Looped only " << nLoopedPlaneIDs
        << " plane IDs, while we expected " << nPlanes << " iterations!";
      ++nErrors;
    } // if
    unsigned int nLoopedPlanes = 0;
    for (geo::PlaneGeo const& Plane: geom->IteratePlanes()) {
      if (nLoopedPlanes >= nPlanes) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "After all " << nLoopedPlanes
          << " planes, iterator has not reached the end";
        ++nErrors;
        break;
      }
      ++nLoopedPlanes;
    }
    if (nLoopedPlanes < nPlanes) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "Looped only " << nLoopedPlanes
        << " planes, while we expected " << nPlanes << " iterations!";
      ++nErrors;
    } // if
    
    if (runningWID) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "Wire ID still valid (" << runningWID
        << ") after incrementing from the last one.";
      ++nErrors;
    }
    
    if (iWireID) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "wire iterator thinks it's still at " << *iWireID
        << ", but we are already finished";
      ++nErrors;
    }
    try {
      geo::WireGeo const& Wire = *iWire;
      LOG_ERROR("GeometryIteratorLoopTest")
        << "wire iterator thinks it's still at " << iWire.ID()
        << ", but we are already finished";
      ++nErrors;
    }
    catch (cet::exception const&) {
      LOG_DEBUG("GeometryIteratorLoopTest") << "exception caught"
        " while dereferencing an iterator to a past-the-end wire.\n";
    }
    
    // test if we can loop all wires with the iterators (via iterator box)
    LOG_DEBUG("GeometryIteratorsDump")
      << "Looping though " << nWires << " wires";
    unsigned int nLoopedWireIDs = 0;
    for (geo::WireID const& wID: geom->IterateWireIDs()) {
      LOG_TRACE("GeometryIteratorsDump") << wID;
      if (nLoopedWireIDs >= nWires) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "After all " << nLoopedWireIDs
          << " wire IDs, iterator has not reached the end ("
          << *(geom->end_wire_id()) << ") but it's still at " << wID;
        ++nErrors;
        break;
      }
      ++nLoopedWireIDs;
    }
    if (nLoopedWireIDs < nWires) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "Looped only " << nLoopedWireIDs
        << " wire IDs, while we expected " << nWires << " iterations!";
      ++nErrors;
    } // if
    unsigned int nLoopedWires = 0;
    for (geo::WireGeo const& Wire: geom->IterateWires()) {
      if (nLoopedWires >= nWires) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "After all " << nLoopedWires
          << " wires, iterator has not reached the end";
        ++nErrors;
        break;
      }
      ++nLoopedWires;
    }
    if (nLoopedWires < nWires) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "Looped only " << nLoopedWires
        << " wires, while we expected " << nWires << " iterations!";
      ++nErrors;
    } // if
    
    
    if (runningSID) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "TPC set ID still valid (" << runningSID
        << ") after incrementing from the last one.";
      ++nErrors;
    }
    
    if (iTPCsetID) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "TPC set iterator thinks it's still at " << *iTPCsetID
        << ", but we are already finished";
      ++nErrors;
    }
    
    // test if we can loop all TPC sets with the iterators (via iterator box)
    LOG_DEBUG("GeometryIteratorsDump")
      << "Looping though " << cumTPCsets << " TPC sets";
    unsigned int nLoopedTPCsetIDs = 0;
    for (readout::TPCsetID const& sID: geom->IterateTPCsetIDs()) {
      LOG_TRACE("GeometryIteratorsDump") << sID;
      if (nLoopedTPCsetIDs >= cumTPCsets) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "After all " << nLoopedTPCsetIDs
          << " TPC set IDs, iterator has not reached the end ("
          << *(geom->end_TPCset_id()) << ") but it's still at " << sID;
        ++nErrors;
        break;
      }
      ++nLoopedTPCsetIDs;
    }
    if (nLoopedTPCsetIDs < cumTPCsets) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "Looped only " << nLoopedTPCsetIDs
        << " TPC set IDs, while we expected " << cumTPCsets << " iterations!";
      ++nErrors;
    } // if
    
    
    if (runningRID) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "readout plane ID still valid (" << runningRID
        << ") after incrementing from the last one.";
      ++nErrors;
    }
    
    if (iROPID) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "readout plane iterator thinks it's still at " << *iROPID
        << ", but we are already finished";
      ++nErrors;
    }
    
    // test if we can loop all planes with the iterators (via iterator box)
    LOG_DEBUG("GeometryIteratorsDump")
      << "Looping though " << cumROPs << " readout planes";
    unsigned int nLoopedROPIDs = 0;
    for (readout::ROPID const& rID: geom->IterateROPIDs()) {
      LOG_TRACE("GeometryIteratorsDump") << rID;
      if (nLoopedROPIDs >= cumROPs) {
        LOG_ERROR("GeometryIteratorLoopTest")
          << "After all " << nLoopedROPIDs
          << " readout planes, ID iterator has not reached the end ("
          << *(geom->end_ROP_id()) << ") but it's still at " << rID;
        ++nErrors;
        break;
      }
      ++nLoopedROPIDs;
    }
    if (nLoopedROPIDs < cumROPs) {
      LOG_ERROR("GeometryIteratorLoopTest")
        << "Looped only " << nLoopedROPIDs
        << " readout plane IDs, while we expected " << cumROPs
        << " iterations!";
      ++nErrors;
    } // if
    
    return nErrors;
  } // GeometryIteratorLoopTestAlg::Run()
  
  //----------------------------------------------------------------------------

#pragma GCC diagnostic pop

} // namespace geo
