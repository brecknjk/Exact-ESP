#ifndef OSM_GRAPH_TOOLS_PROCESSORS_H
#define OSM_GRAPH_TOOLS_PROCESSORS_H

#include <unordered_map>
#include <sstream>
#include <chrono>

#include "../utils/types.hpp"

#ifdef USE_OSMPBF
	#include "osmpbf/pbistream.h"
	#include <osmpbf/iway.h>
	#include <osmpbf/inode.h>
	#include <osmpbf/filter.h>
	#include <osmpbf/primitiveblockinputadaptor.h>

	void gatherEdgesAndNodes(std::vector<std::string> & inputFileNames, std::vector<Coordinates> & coords, std::vector<Edge> & edges) {
		osmpbf::PbiStream inFile = osmpbf::PbiStream(inputFileNames);
		osmpbf::PrimitiveBlockInputAdaptor pbi;
		std::unordered_map<int64_t, uint32_t> keepNodes;

		// gather edges
		inFile.dataSeek(0);
		while (inFile.parseNextBlock(pbi)) {
			if (!pbi.isNull() && pbi.waysSize()) {
				uint32_t naturalTagId = pbi.findString("natural");
				if (naturalTagId == 0) continue;
				for (osmpbf::IWayStream way = pbi.getWayStream(); !way.isNull(); way.next()) {
					if ("coastline" != way.valueByKeyId(naturalTagId)) continue;

					uint32_t lastLocalID = UINT32_MAX;
					for(osmpbf::IWayStream::RefIterator refIt(way.refBegin()), refEnd(way.refEnd()); refIt != refEnd; ++refIt) {
						if(keepNodes.count(*refIt) == 0) keepNodes[*refIt] = keepNodes.size();
						if(lastLocalID != UINT32_MAX) {
							edges.push_back(Edge(lastLocalID, keepNodes[*refIt]));
						}
						lastLocalID = keepNodes[*refIt];
					}
				}
			}
		}
		// gather nodes
		uint32_t foundNodes = 0;
		coords.resize(keepNodes.size());
		inFile.dataSeek(0);
		while (inFile.parseNextBlock(pbi)) {
			if (!pbi.isNull() && pbi.nodesSize()) {
				for (osmpbf::INodeStream node = pbi.getNodeStream(); !node.isNull(); node.next()) {
					if (keepNodes.count(node.id()) > 0){
						coords[keepNodes[node.id()]] = Coordinates(node.latd(), node.lond());
					}
				}
			}
		}
	}
#else
	void gatherEdgesAndNodes(std::vector<std::string> & inputFileNames, std::vector<Coordinates> & coords, std::vector<Edge> & edges) {
		throw std::runtime_error("No OSM PBF support compiled in, please compile with -DUSE_OSMPBF=ON");
	}
#endif

#endif
