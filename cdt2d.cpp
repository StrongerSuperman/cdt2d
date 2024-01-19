#include <cassert>
#include <iostream>

#include "predicates.hpp"
#include "cdt2d.hpp"


#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

namespace CDT
{
	Triangulation::Triangulation()
	{
		goodAngleCosSquare = cos(minAngle * PI / 180.0);
		goodAngleCosSquare *= goodAngleCosSquare;
		if (goodAngleCosSquare == 1.0)
		{
			offConstant = 0.0;
		}
		else
		{
			offConstant = 0.475 * sqrt((1.0 + goodAngleCosSquare) / (1.0 - goodAngleCosSquare));
		}
	}

	void Triangulation::InsertVertices(const std::vector<Vertex>& vertexList)
	{
		this->vertices = vertexList;
	}

	void Triangulation::InsertEdges(const std::vector<Edge>& edgeList, const std::vector<std::vector<Edge>>& holesEdgeList)
	{
		this->edges = edgeList;
		this->holesEdges = holesEdgeList;
	}

	void Triangulation::Perform(bool delaunayRefine)
	{
		std::cout << "------------------Start--------------------" << std::endl;
		size_t verticesSize = vertices.size();
		assert(vertices.size() > 2);

		// init triangulation by first two vertices and a ghost vertex
		IdxType firstVtx = 0;
		IdxType secondVtx = 1;
		addTriangle(GHOST_VTX_IDX, firstVtx, secondVtx);
		addTriangle(GHOST_VTX_IDX, secondVtx, firstVtx);
		mostRecentAddedTri = getTriangleIdx(Edge(firstVtx, secondVtx));

		// insert vertex
		for (IdxType v = secondVtx + 1; v < verticesSize; v++)
		{
			insertVertex(v);
		}

		// insert edge
		for (const Edge& edge : this->edges)
		{
			insertEdge(edge);
		}
		// NOTE: the holeEdge's two vertices must be in CCW
		//for (const auto& holeEdges : this->holesEdges)
		//{
		//	for (const Edge& holeEdge : holeEdges)
		//	{
		//		insertEdge(holeEdge);
		//	}
		//}

		removeOuterAndHoles();

		remapTriangulation();

		if (delaunayRefine)
		{
			delaunayRefinement();
			remapTriangulation();
		}

		std::cout << "total triangles num: " << triangles.size() << std::endl;
		std::cout << "-------------------End---------------------" << std::endl;
		std::cout << std::endl;
	}

	void Triangulation::insertVertex(IdxType v)
	{
		/* The implementation is based on 'Bowyer-Watson' incremental insertion algorithm
		 * In expectation, fewer than six new triangles created for one vertex insertion
		 * Point location expected time complexity: O(NLogN)
		 * Vertex insertion expected time complexity: O(N)
		 * 
		 * TODO: use 'BRIO' to guarantee the expected running time.
		 * A biased randomized insertion order (BRIO) is a permutation of the vertices that has strong spatial
		 * locality but retains enough randomness to obtain an expected running time in O(n2). Their experiments
		 * show that a BRIO greatly improves the efficiency of the memory hierarchy��especially virtual memory.
		 * 
		 * There two 'BRIO' methods:
		 * 1.insert the vertices in the order they are encountered on a space-filling curve such as Hilbert curve
		 * or a z-order curve.
		 * 2.store the vertices in an octree or k-d tree, refined so each leaf node contains only a few vertices,
		 * then the vertices of a traversal of the tree
		 */
		TriangulationLocation loc;
		locatePoint(v, loc);
		if (loc.type == TriangulationLocationType::TL_INSIDE)
		{
			insertVtxInsideTriangulation(v, loc.tri);
		}
		else if(loc.type == TriangulationLocationType::TL_OUTSIDE)
		{
			insertVtxOutsideTriangulation(v, loc);
		}
		else
		{
			assert(false);
		}
	}

	void Triangulation::locatePoint(IdxType v, TriangulationLocation& loc)
	{
		/* Walking from a specific triangle to the triangle that contains v
		 * 
		 * Walking is fast in practice if it follows two guides:
		 * 1.the vertices should be inserted in an order that has much spatial locality
		 * 2.each walk should begin at the most recently created triangle
		 * 
		 * The typical walk visits small constant number of triangle
		 */

		assert(edgeTriIdxTable.size() > 3);   // at least three ghost triangles exist

		PntLineLocationType plLocType;
		TriUSet visited;
		std::stack<const Triangle*> candidates;
		candidates.push(getTriangle(mostRecentAddedTri));
		mostRecentAddedTri = DUMMY_TRI_IDX;
		const Triangle* neighborTri = nullptr;
		IdxType vStart;
		IdxType vEnd;
		while (!candidates.empty())
		{
			const Triangle* curTri = candidates.top();
			candidates.pop();
			assert(curTri != nullptr);

			// NOTE: the triangulation of a point set is a convex hull
			if (isTriangleGhost(*curTri))
			{
				vStart = curTri->GetVtxCCW(GHOST_VTX_IDX);
				vEnd = curTri->GetVtxCW(GHOST_VTX_IDX);
				plLocType = locatePntLine(v, vStart, vEnd);
				if (plLocType == PntLineLocationType::PL_LEFT ||
					plLocType == PntLineLocationType::PL_ON_LINE)
				{
					loc.type = TriangulationLocationType::TL_OUTSIDE;
					break;
				}
				neighborTri = getTriangle(Edge(vEnd, vStart));
				assert(neighborTri != nullptr);
				candidates.push(neighborTri);
			}
			else
			{
				bool foundInTriangle = true;
				bool foundInBoundary = false;
				// TODO: use stochastic offset to randomize which edge we check first
				for (IdxType i = 0; i < 3; ++i)
				{
					vStart = curTri->V(i % 3);
					vEnd = curTri->V((i + 1) % 3);
					neighborTri = getTriangle(Edge(vEnd, vStart));
					assert(neighborTri != nullptr);
					plLocType = locatePntLine(v, vStart, vEnd);
					if (plLocType == PntLineLocationType::PL_RIGHT &&
						visited.insert(*neighborTri).second)
					{
						foundInTriangle = false;
						candidates.push(neighborTri);
						break;
					}
					else if (plLocType == PntLineLocationType::PL_ON_LINE &&
						isTriangleGhost(*neighborTri))
					{
						foundInBoundary = true;
						loc.type = TriangulationLocationType::TL_OUTSIDE;
						break;
					}
				}

				if (foundInBoundary)
				{
					break;
				}
				if (foundInTriangle)
				{
					loc.type = TriangulationLocationType::TL_INSIDE;
					loc.tri = curTri;
					break;
				}
			}
		}

		if (loc.type == TriangulationLocationType::TL_OUTSIDE)
		{
			// collect all outer edges which right half plane enclose the vertex
			for (const Edge& edge : ghostTriEdges)
			{
				if (isEdgeGhost(edge))
					continue;
				plLocType = locatePntLine(v, edge.V1(), edge.V2());
				if (plLocType == PntLineLocationType::PL_LEFT)
				{
					loc.inHPEdges.push_back(edge);
				}
				else if (plLocType == PntLineLocationType::PL_ON_LINE)
				{
					loc.onHPEdges.push_back(edge);
					// check vertex v whether in middle of the edge
					if (isPntInMidOfLine(v, edge.V1(), edge.V2()))
					{
						loc.boundaryEdge = edge;
						loc.onBoundaryEdge = true;
						break;
					}
				}
			}
		}
	}

	void Triangulation::insertVtxInsideTriangulation(IdxType v, const Triangle* tri)
	{
		PntTriLocationType locTri = locatePntTriangle(v, tri->V1(), tri->V2(), tri->V3());
		assert(locTri != PntTriLocationType::PT_OUTSIDE);
		if (locTri > PntTriLocationType::PT_INSIDE)  /*on edge*/
		{
			size_t edgeIdxInTri = static_cast<size_t>(locTri - PntTriLocationType::PT_ON_EDGE_1);
			size_t v1IdxInTri = edgeIdxInTri;
			size_t v2IdxInTri = (v1IdxInTri + 1) % 3;
			insertVtxOnEdge(v, tri->V(v1IdxInTri), tri->V(v2IdxInTri));
		}
		else  /*inside triangle*/
		{
			insertVtxInTriangle(v, tri);
		}
	}

	void Triangulation::insertVtxOnEdge(IdxType v, IdxType edgV1, IdxType edgV2)
	{
		Edge edge(edgV1, edgV2);
		Edge edgeRev(edgV2, edgV1);
		IdxType edgTriV = adjacent(edgV1, edgV2);
		IdxType edgRevTriV = adjacent(edgV2, edgV1);
		deleteTriangle(edgV1, edgV2, edgTriV);
		deleteTriangle(edgV2, edgV1, edgRevTriV);
		digCavity(v, edgV2, edgTriV);
		digCavity(v, edgTriV, edgV1);
		digCavity(v, edgV1, edgRevTriV);
		digCavity(v, edgRevTriV, edgV2);
	}

	void Triangulation::insertVtxInTriangle(IdxType v, const Triangle* tri)
	{
		IdxType v1 = tri->V1();
		IdxType v2 = tri->V2();
		IdxType v3 = tri->V3();
		deleteTriangle(v1, v2, v3);
		digCavity(v, v1, v2);
		digCavity(v, v2, v3);
		digCavity(v, v3, v1);
	}

	void Triangulation::insertVtxOutsideTriangulation(IdxType v, const TriangulationLocation& loc)
	{
		/* NOTE: the vertex v is on the left of these edges*/

		if (loc.onBoundaryEdge)
		{
			const Edge& boundaryEdge = loc.boundaryEdge;
			assert(isVtxValid(boundaryEdge.V1()) && isVtxValid(boundaryEdge.V2()));
			insertVtxOnEdge(v, boundaryEdge.V1(), boundaryEdge.V2());
		}
		else
		{
			/* NOTE: Since delaunay triangulation of a point set is a convex, so we can conclude that
			the edge in inHPEdges or onHPEdges connect each other(except start and end edge of the edge chain)*/

			assert(!(loc.inHPEdges.empty() && loc.onHPEdges.empty()));

			auto findStartEndVertices = [](const std::vector<Edge>& edges)
			{
				std::unordered_map<IdxType, int> vtx2CountMap;
				for (const Edge& edge : edges)
				{
					std::vector<IdxType> vertices({ edge.V1(), edge.V2() });
					for (IdxType vertex : vertices)
					{
						if (vtx2CountMap.find(vertex) == vtx2CountMap.end())
						{
							vtx2CountMap[vertex] = 1;
						}
						else
						{
							vtx2CountMap[vertex] += 1;
						}
					}
				}
				std::vector<IdxType> countOneVertices;
				for (auto& vtx2CountIter : vtx2CountMap)
				{
					if (vtx2CountIter.second > 1)
						continue;
					countOneVertices.push_back(vtx2CountIter.first);
				}
				assert(countOneVertices.size() == 2);
				IdxType v1 = countOneVertices[0];
				IdxType v2 = countOneVertices[1];
				return std::make_pair(v1, v2);
			};

			if (!loc.inHPEdges.empty())
			{
				auto& edgeChains = loc.inHPEdges;
				// delete old ghost triangles
				for (const Edge& edge : edgeChains)
				{
					deleteTriangle(edge.V1(), edge.V2(), GHOST_VTX_IDX);
				}
				// add new solid triangles and delete some old solid triangles
				for (const Edge& edge : edgeChains)
				{
					digCavity(v, edge.V1(), edge.V2());
				}
				// add new ghost triangles (TODO: use a elegant way to handle this)
				std::pair<IdxType, IdxType> startEndVertices = findStartEndVertices(edgeChains);
				IdxType v1 = startEndVertices.first;
				IdxType v2 = startEndVertices.second;
				PntLineLocationType plLocType = locatePntLine(v, v1, v2);
				if (plLocType == PntLineLocationType::PL_LEFT)
				{
					addTriangle(GHOST_VTX_IDX, v1, v);
					addTriangle(GHOST_VTX_IDX, v, v2);
				}
				else if (plLocType == PntLineLocationType::PL_RIGHT)
				{
					addTriangle(GHOST_VTX_IDX, v2, v);
					addTriangle(GHOST_VTX_IDX, v, v1);
				}
				else
				{
					assert(false);
				}
			}
			else
			{
				/* the following branch only occur when all inserted vertices are colliner with a line*/

				// remove the same edge only with different direction
				std::vector<Edge> edgeChains;
				for (const Edge& edgeOrg : loc.onHPEdges)
				{
					bool skip = false;
					for (const Edge& edge : edgeChains)
					{
						if (edgeOrg == edge || Edge(edgeOrg.V2(), edgeOrg.V1()) == edge)
						{
							skip = true;
						}
					}
					if (!skip)
					{
						edgeChains.push_back(edgeOrg);
					}
				}
				// find the closet vertex of these edge to vertex v
				std::pair<IdxType, IdxType> StartEndVertices = findStartEndVertices(edgeChains);
				IdxType v1 = StartEndVertices.first;
				IdxType v2 = StartEndVertices.second;
				Vertex vtx1 = vertices[v1];
				Vertex vtx2 = vertices[v2];
				Vertex vtx = vertices[v];
				PrecisionType distX1 = vtx1.X() - vtx.X();
				PrecisionType distY1 = vtx1.Y() - vtx.Y();
				PrecisionType distX2 = vtx2.X() - vtx.X();
				PrecisionType distY2 = vtx2.Y() - vtx.Y();
				PrecisionType dist1Square = distX1 * distX1 + distY1 * distY1;
				PrecisionType dist2Square = distX2 * distX2 + distY2 * distY2;
				IdxType closetV = dist1Square > dist2Square ? v2 : v1;
				// find the edge that the closet point on
				const Edge* closetEdge;
				for (const Edge& edge : edges)
				{
					if (edge.V1() == closetV || edge.V2() == closetV)
					{
						closetEdge = &edge;
						break;
					}
				}
				// triangulate
				addTriangle(GHOST_VTX_IDX, closetV, v);
				addTriangle(GHOST_VTX_IDX, v, closetV);
			}
		}
	}

	void Triangulation::digCavity(IdxType v, IdxType v1, IdxType v2)
	{
		if (isTriangleGhost(v, v1, v2))
		{
			addTriangle(v, v1, v2);
			return;
		}

		IdxType vOpo = adjacent(v2, v1);
		if (!isVtxValid(vOpo))  // do nothing if the triangle has already been deleted
			return;

		/* depth-first search to find all triangles that are no longer Delaunay*/
		if (isVtxGhost(vOpo) || isTriangleDelaunay(v, v1, v2, vOpo))
		{
			addTriangle(v, v1, v2);
		}
		else
		{
			deleteTriangle(v2, v1, vOpo);
			digCavity(v, v1, vOpo);
			digCavity(v, vOpo, v2);
		}
	}

	void Triangulation::insertEdge(const Edge& edge)
	{
		/* NOTE:
		 * 1. The edge vertex must be already inserted into the triangulaton
		 * 2. The handle for INTERSECTION case between two edge is not support now
		*/

		// check if the insert edge is already the edge of a triangle in this triangulation
		if (isEdgeTriExist(edge))
		{
			addConstrainEdge(edge);
			return;
		}

		Edge insertEdge = edge;
		std::vector<Edge> remainingEdges;
		remainingEdges.push_back(insertEdge);
		while (!remainingEdges.empty())
		{
			insertEdge = remainingEdges.back();
			remainingEdges.pop_back();
			insertEdgeIter(insertEdge, remainingEdges);
		}
	}

	void Triangulation::insertEdgeIter(const Edge& edge, std::vector<Edge>& remainingEdges)
	{
		assert(!isEdgeGhost(edge));

		IdxType vEdgeStart = edge.V1();
		IdxType vEdgeEnd = edge.V2();
		if (vEdgeStart == vEdgeEnd)   // skip edge that connect a vertex to itself
			return;

		// find some candidate triangles in edge start vertex
		const Triangle* curTri = nullptr;
		IdxType v = vEdgeStart;
		IdxType vLeft;
		IdxType vRight;
		std::tie(vRight, vLeft) = adjacentOne(v);
		curTri = getTriangle(Edge(vLeft, vEdgeStart));
		std::vector<const Triangle*> candidateTris(1, curTri);
		IdxType vLeftOpo = getOppositeVtx(*curTri, vLeft);
		IdxType vCurTriLeft = vLeft;
		while (vLeftOpo != vLeft)
		{
			vCurTriLeft = curTri->GetVtxCW(vEdgeStart);
			curTri = getOppositeTriangle(*curTri, vCurTriLeft);
			if (curTri == nullptr)
				break;
			candidateTris.push_back(curTri);
			vCurTriLeft = curTri->GetVtxCW(vEdgeStart);
			vLeftOpo = getOppositeVtx(*curTri, vCurTriLeft);
		}
		if (curTri != nullptr)
			curTri = getOppositeTriangle(*curTri, vCurTriLeft);
		if (curTri != nullptr)
			candidateTris.push_back(curTri);

		// find intersected triangle in  edge start vertex
		PntLineLocationType plLocType;
		curTri = nullptr;   // reset curTri
		for (const Triangle* canTri : candidateTris)
		{
			/* NOTE: the intersect triangle can be a triangle that
			its two edge vEdgeStart-vLeft and vEdgeStart-vRight collinear with vEdgeStart-vEdgeEnd*/

			// TODO: does edge really never start from a ghost triangle?
			if (isTriangleGhost(*canTri))
				continue;
			vLeft = canTri->GetVtxCW(vEdgeStart);
			vRight = canTri->GetVtxCCW(vEdgeStart);
			plLocType = locatePntLine(vLeft, vEdgeStart, vEdgeEnd);
			if (plLocType == PntLineLocationType::PL_LEFT)
			{
				plLocType = locatePntLine(vRight, vEdgeStart, vEdgeEnd);
				if (plLocType == PntLineLocationType::PL_RIGHT)
				{
					// intersect triangle found
					curTri = canTri;
					break;
				}
				else if (plLocType == PntLineLocationType::PL_ON_LINE)
				{
					/* if one of the triangle vertices is on the edge, move edge start*/

					// check whether vEdgeStart-vLeft is on the same side with vEdgeStart-vEdgeEnd
					if (!isPntInMidOfLine(vEdgeStart, vRight, vEdgeEnd))
					{
						Edge edgePart(vEdgeStart, vRight);
						addConstrainEdge(edgePart);
						Edge edgeRemain(vRight, vEdgeEnd);
						remainingEdges.push_back(edgeRemain);
						return;
					}
				}
			}
			else if (plLocType == PntLineLocationType::PL_ON_LINE)
			{
				/* if one of the triangle vertices is on the edge, move edge start*/

				// check whether vEdgeStart-vLeft is on the same side with vEdgeStart-vEdgeEnd
				if (!isPntInMidOfLine(vEdgeStart, vLeft, vEdgeEnd))
				{
					Edge edgePart(vEdgeStart, vLeft);
					addConstrainEdge(edgePart);
					Edge edgeRemain(vLeft, vEdgeEnd);
					remainingEdges.push_back(edgeRemain);
					return;
				}
			}
		}

		assert(curTri != nullptr);

		// walk from edge start vertex to end to find all triangles that intersect the edge
		IdxType vOpo;
		const Triangle* curTriOpo;
		std::vector<IdxType> ptsLeft(1, vLeft);
		std::vector<IdxType> ptsRight(1, vRight);
		std::vector<const Triangle*> isectTris(1, curTri);
		IdxType vCurEnd = vEdgeEnd;
		while (!curTri->IsContainVtx(vCurEnd))
		{
			curTriOpo = getOppositeTriangle(*curTri, v);
			vOpo = getOppositeVtx(*curTri, v);

			assert(isVtxValid(vOpo));
			assert(!isVtxGhost(vOpo));

			isectTris.push_back(curTriOpo);
			curTri = curTriOpo;

			plLocType = locatePntLine(vOpo, vEdgeStart, vEdgeEnd);
			if (plLocType == PntLineLocationType::PL_LEFT)
			{
				ptsLeft.push_back(vOpo);
				v = vLeft;
				vLeft = vOpo;
			}
			else if (plLocType == PntLineLocationType::PL_RIGHT)
			{
				ptsRight.push_back(vOpo);
				v = vRight;
				vRight = vOpo;
			}
			else
			{
				vCurEnd = vOpo;
			}
		}

		// delete triangles that intersect the edge
		for (const Triangle* tri : isectTris)
		{
			deleteTriangle(tri->V1(), tri->V2(), tri->V3());
		}

		// retriangulate left and right pseudo polygons
		triangulatePseudopolygon(vEdgeStart, vCurEnd, std::move(ptsLeft));
		std::reverse(ptsRight.begin(), ptsRight.end());
		triangulatePseudopolygon(vCurEnd, vEdgeStart, std::move(ptsRight));

		// check whether walking reach at the end of the edge
		if (vCurEnd != vEdgeEnd) // encountered point on the edge
		{
			Edge edgePart(vEdgeStart, vCurEnd);
			addConstrainEdge(edgePart);
			Edge edgeRemain(vCurEnd, vEdgeEnd);
			remainingEdges.push_back(edgeRemain);
			return;
		}
		else
		{
			addConstrainEdge(edge);
		}
	}

	void Triangulation::triangulatePseudopolygon(IdxType vStart, IdxType vEnd, std::vector<IdxType>&& midPts)
	{
		assert(!midPts.empty());
		typedef std::tuple<IdxType, IdxType, std::vector<IdxType>> Pseudopolygon;
		std::vector<Pseudopolygon> remaining;
		remaining.emplace_back(vStart, vEnd, std::move(midPts));
		IdxType vCurStart = vStart;
		IdxType vCurEnd = vEnd;
		std::vector<IdxType> curMidPts;
		while (!remaining.empty())
		{
			std::tie(vCurStart, vCurEnd, curMidPts) = remaining.back();
			remaining.pop_back();
			// check if pseudo-polygon is a single triangle
			if (curMidPts.size() == 1)
			{
				addTriangle(vCurStart, vCurEnd, curMidPts.back());
				continue;
			}
			// find Delaunay point
			size_t vDelaunayIdx = findDelaunayPoint(vCurStart, vCurEnd, curMidPts);
			IdxType vDelaunay = curMidPts[vDelaunayIdx];
			// add new triangle
			addTriangle(vCurStart, vCurEnd, vDelaunay);
			// split pseudo-polygon in two parts and triangulate them
			std::vector<IdxType> curMidPts1(curMidPts.begin(), curMidPts.begin() + vDelaunayIdx);
			if (!curMidPts1.empty())
			{
				remaining.emplace_back(vCurStart, vDelaunay, std::move(curMidPts1));
			}
			std::vector<IdxType> curMidPts2(curMidPts.begin() + vDelaunayIdx + 1, curMidPts.end());
			if (!curMidPts2.empty())
			{
				remaining.emplace_back(vDelaunay, vCurEnd, std::move(curMidPts2));
			}
		}
	}

	size_t Triangulation::findDelaunayPoint(IdxType vStart, IdxType vEnd, const std::vector<IdxType>& midPts)
	{
		assert(!midPts.empty());
		size_t vDelaunayPtIdx = 0;
		for (size_t i = 0; i < midPts.size(); i++)
		{
			if (!isTriangleDelaunay(vStart, vEnd, midPts[vDelaunayPtIdx], midPts[i]))
			{
				vDelaunayPtIdx = i;
			}
		}
		return vDelaunayPtIdx;
	}

	void Triangulation::removeOuterAndHoles()
	{
		if (0)
		{
			// only reomve outer(ghost triangles)
			auto iter = edgeTriIdxTable.begin();
			while (iter != edgeTriIdxTable.end())
			{
				if (isTriangleGhost(*getTriangle(iter->second)))
					iter = edgeTriIdxTable.erase(iter);
				else
					iter++;
			}
			return;
		}

		std::vector<LayerDepth> triDepths;
		calculateTriangleDepths(triDepths);
		TriIdxUSet toErase;
		toErase.reserve(triangles.size());
		for (size_t idx = 0; idx < triangles.size(); ++idx)
		{
			if (triDepths[idx] % 2 == 0)
				toErase.insert(static_cast<IdxType>(idx));
		}
		auto iter = edgeTriIdxTable.begin();
		while (iter != edgeTriIdxTable.end())
		{
			if (toErase.count(iter->second))
				iter = edgeTriIdxTable.erase(iter);
			else
				iter++;
		}
	}

	void Triangulation::calculateTriangleDepths(std::vector<LayerDepth>& triDepths)
	{
		/**
		 * Calculate depth of each triangle in constraint triangulation.
		 *
		 * Perform depth peeling from super triangle to outermost boundary,
		 * then to next boundary and so on until all triangles are traversed.@n
		 * For example depth is:
		 * 0 for triangles outside outermost boundary
		 * 1 for triangles inside boundary but outside hole
		 * 2 for triangles in hole
		 * 3 for triangles in island and so on...
		 */
		triDepths = std::vector<LayerDepth>(triangles.size(), DEFAULT_LAYER_DEPTH);
		IdxType vCCW;
		IdxType vCW;
		std::tie(vCCW, vCW) = adjacentOne(GHOST_VTX_IDX);
		IdxType triIdx = getTriangleIdx(Edge(vCCW, vCW));
		TriIdxStack seeds;
		seeds.push(triIdx);
		LayerDepth layerDepth = 0;
		LayerDepth deepestSeedDepth = 0;
		LayerTriIdxesUMap seedsByDepth;
		do
		{
			const TriLayerUMap& newSeeds = peelLayer(seeds, layerDepth, triDepths);
			seedsByDepth.erase(layerDepth);
			for (auto& triLayerPair : newSeeds)
			{
				deepestSeedDepth = deepestSeedDepth > triLayerPair.second ? deepestSeedDepth : triLayerPair.second;
				seedsByDepth[triLayerPair.second].insert(triLayerPair.first);
			}
			const TriIdxUSet& nextLayerSeeds = seedsByDepth[layerDepth + 1];
			seeds = TriIdxStack(TriIdxDeque(nextLayerSeeds.begin(), nextLayerSeeds.end()));
			++layerDepth;
		} while (!seeds.empty() || deepestSeedDepth > layerDepth);
	}

	TriLayerUMap Triangulation::peelLayer(
		TriIdxStack seeds, const LayerDepth layerDepth, std::vector<LayerDepth>& triDepths)
	{
		/**
		 * Depth-peel a layer in triangulation, used when calculating triangle depths
		 *
		 * It takes starting seed triangles, traverses neighboring triangles, and
		 * assigns given layer depth to the traversed triangles. Traversal is
		 * blocked by constraint edges. Triangles behind constraint edges are
		 * recorded as seeds of next layer and returned from the function.
		 */
		TriLayerUMap behindBoundary;
		while (!seeds.empty())
		{
			const IdxType triIdx = seeds.top();
			seeds.pop();
			triDepths[triIdx] = layerDepth;
			behindBoundary.erase(triIdx);
			const Triangle& tri = triangles[triIdx];
			for (size_t i = 0; i < 3; ++i)
			{
				IdxType v = tri.V(i);
				IdxType triOpoIdx = getOppositeTriangleIdx(tri, v);
				if (triOpoIdx == DUMMY_TRI_IDX || triDepths[triOpoIdx] <= layerDepth)
					continue;
				Edge opEdge(tri.GetVtxCCW(v), tri.GetVtxCW(v));
				if (isEdgeConstrained(opEdge))
				{
					const LayerDepth triDepth = layerDepth + 1;
					behindBoundary[triOpoIdx] = triDepth;
					continue;
				}
				seeds.push(triOpoIdx);
			}
		}
		return behindBoundary;
	}

	void Triangulation::remapTriangulation()
	{
		std::vector<Triangle> triangleRemap;
		triangleRemap.reserve((triangles.size()));
		std::unordered_map<Triangle, IdxType, HashTriangle> triIdxMap;
		int i = 0;
		for (auto& edgeTriIdxPair : edgeTriIdxTable)
		{
			const Edge& edge = edgeTriIdxPair.first;
			IdxType triIdx = edgeTriIdxPair.second;
			const Triangle& tri = triangles[triIdx];
			auto iter = triIdxMap.find(tri);
			if (iter == triIdxMap.end())
			{
				triIdxMap.emplace(tri, i);
				edgeTriIdxTable[edge] = i;
				triangleRemap.push_back(tri);
				i++;
			}
			else
			{
				edgeTriIdxTable[edge] = iter->second;
			}
		}
		triangleRemap.erase(triangleRemap.begin() + i, triangleRemap.end());
		std::swap(triangles, triangleRemap);
	}

	void Triangulation::delaunayRefinement()
	{
		/* Delaunay refinment based on Ruppert's algorithm.
         
         * The algorithm interleaves segment splitting with triangle splitting,
         * and encroached segments are given priority over skinny triangles.
		 * 
		 * Step1: Segment Split
		 * Desc:  Find all 'encroached' segments to split(segment here is the constrained edge in CDT).
		 *        A subsegment is said to be encroached if a vertex other than its endpoints lies on or
		 *        inside its diametral circle, and the encroaching vertex is visible from the interior
		 *        of the segments(Visibility is obstructed only by other segments).
		 * 
		 * Step2: Triangle Split
		 * Desc:  Find bad triangles and insert point in its circumcenter.
		 * 
		 * NOTE:  The input is CDT, the outer triangles and holes have already been removed.
		 */
        
        insertSteinerCnt = 0;
        failedSplitBadTris = 0;
        totalSplitTris = 0;
        totalSplitEdges = 0;
        EdgeStack().swap(encroachedEdges);
        badTriangles.clear();

		// collect all encroached edges
		EdgeUSet checkedEdges;
		for (const auto& edgeTriPair : edgeTriIdxTable)
		{
			const Edge& edge = edgeTriPair.first;
			checkedEdges.insert(edge);
			Edge edgeRev(edge.V2(), edge.V1());
			if (checkedEdges.find(edgeRev) != checkedEdges.end())
				continue;
			if (isEdgeEncroached(edge))
			{
				encroachedEdges.push(edge);
			}
		}
		checkedEdges.clear();
        
        // fix encroached edges without noting bad triangles
        while(!encroachedEdges.empty() && insertSteinerCnt < maxSteinerCnt)
        {
            Edge edge = encroachedEdges.top();
            encroachedEdges.pop();
            
            splitEdge(edge);
        }

		// collect all bad triangles
		for (size_t i = 0; i < triangles.size(); i++)
		{
			const Triangle& tri = triangles[i];
			TriangleQuality triQuality;
			calTriangleQuality(tri, triQuality);
			// Check whether the angle is smaller than permitted
			if (triQuality.minAngleCosSquare > goodAngleCosSquare)
			{
				badTriangles.push_back(i);
			}
		}

        assert(encroachedEdges.empty());
		while (!badTriangles.empty() && insertSteinerCnt < maxSteinerCnt)
		{
            IdxType badTriIdx = badTriangles.front();
            badTriangles.pop_front();
            
            splitTriangle(badTriIdx);
            
            if (!encroachedEdges.empty())
            {
                // put bad triangle back in queue for another try later
                badTriangles.push_back(badTriIdx);
                
                // fix any encroached subsegments that resulted
                while(!encroachedEdges.empty() && insertSteinerCnt < maxSteinerCnt)
                {
                    Edge edge = encroachedEdges.top();
                    encroachedEdges.pop();
                    
                    splitEdge(edge);
                }
            }
		}
        
		std::cout << "total inserted steiner points num: " << insertSteinerCnt << std::endl;
		std::cout << "total splitted encroached edges num: " << totalSplitEdges << std::endl;
		std::cout << "total splitted bad triangles num: " << totalSplitTris << std::endl;
		std::cout << "faild split bad triangles num: " << failedSplitBadTris << std::endl;
		std::cout << std::endl;
	}

	void Triangulation::splitEdge(const Edge& edge)
	{
		const Triangle* tri = getTriangle(edge);
		IdxType edgeV1 = edge.V1();
		IdxType edgeV2 = edge.V2();
		const Vertex& vtx1 = vertices[edgeV1];
		const Vertex& vtx2 = vertices[edgeV2];
		/* To decide where to split a segment, we need to know if the
		 * segment shares an endpoint with an adjacent segment.
		 *
		 * The following strategy is derived from 'Triangle-Library'
		 *
		 * The concern is that, if we simply split every encroached
		 * segment in its center, two adjacent segments with a small
		 * angle between them might lead to an infinite loop; each
		 * vertex added to split one segment will encroach upon the
		 * other segment, which must then be split with a vertex that
		 * will encroach upon the first segment, and so on forever.
		 *
		 * To avoid this, imagine a set of concentric circles, whose
		 * radii are powers of two, about each segment endpoint.
		 * These concentric circles determine where the segment is
		 * split.  (If both endpoints are shared with adjacent
		 * segments, split the segment in the middle, and apply the
		 * concentric circles for later splittings.)
		 */
		PrecisionType split = 0.5;
		IdxType v3 = tri->GetVtxCCW(edgeV2);
		bool acuteorg = isEdgeTriExist(Edge(v3, edgeV2));
		bool acutedest = isEdgeTriExist(Edge(edgeV1, v3));
		if (acuteorg || acutedest)
		{
			PrecisionType segmentlength = sqrt((vtx1.X() - vtx2.X()) * (vtx1.X() - vtx2.X()) +
				(vtx1.Y() - vtx2.Y()) * (vtx1.Y() - vtx2.Y()));
			PrecisionType nearestpoweroftwo = 1.0;
			while (segmentlength > 3.0 * nearestpoweroftwo) {
				nearestpoweroftwo *= 2.0;
			}
			while (segmentlength < 1.5 * nearestpoweroftwo) {
				nearestpoweroftwo *= 0.5;
			}
			// where do we split the segment?
			split = nearestpoweroftwo / segmentlength;
			if (acutedest)
			{
				split = 1.0 - split;
			}
		}

		// create new vertex
		PrecisionType newVtxX = vtx1.X() + (vtx2.X() - vtx1.X()) * split;
		PrecisionType newVtxY = vtx1.Y() + (vtx2.Y() - vtx1.Y()) * split;
		Vertex newVtx(newVtxX, newVtxY);
        IdxType newVtxIdx = vertices.size();
        vertices.push_back(newVtx);

		// try to insert the new vertex on edge
		try
		{
			insertVtxOnEdgeCDT(newVtxIdx, edgeV1, edgeV2);
            
            // check whether the splitted sub-edge encroached
            Edge splitEdge1(edgeV1, newVtxIdx);
            Edge splitEdge2(newVtxIdx, edgeV2);
            if (isEdgeEncroached(splitEdge1))
            {
                encroachedEdges.push(splitEdge1);
            }
            if (isEdgeEncroached(splitEdge2))
            {
                encroachedEdges.push(splitEdge2);
            }

            totalSplitEdges++;
            insertSteinerCnt++;
		}
		catch (...)
		{
            /* something unexpected happended*/
            vertices.pop_back();
		}
	}

	void Triangulation::splitTriangle(IdxType badTriIdx)
	{
        const Triangle& badTri = triangles[badTriIdx];
        
		// add it to vertices for point location operation below
		Vertex centerVtx = findCircumcenter(badTri);
		//Vertex centerVtx = findCentriod(badTri);
        IdxType newVtxIdx = vertices.size();
        vertices.push_back(centerVtx);

		// the new vertex must be in the interior
		TriangulationLocation loc;
		locatePointWithGuide(newVtxIdx, &badTri, loc);
        assert (loc.type == TriangulationLocationType::TL_INSIDE);

		// check steiner vertex too close or duplicate with exist vertex
		const Triangle* tri = loc.tri;
		std::vector<IdxType> triVertices({ tri->V1(), tri->V2(), tri->V3()});
		for (IdxType triV : triVertices)
		{
			const Vertex& vtx = vertices[triV];
			double dxSquare = (vtx.X() - centerVtx.X()) * (vtx.X() - centerVtx.X());
			double dySquare = (vtx.Y() - centerVtx.Y())* (vtx.Y() - centerVtx.Y());
			if (dxSquare + dySquare < CDT_ZERO * CDT_ZERO)
			{
				failedSplitBadTris++;
				vertices.pop_back();
				return;
			}
		}

		// try to insert the new vertex
		try
		{
            VertexInsertType insertResult = insertVtxCDT(newVtxIdx, tri);
            if(insertResult == VI_SUCCESS)
            {
                insertSteinerCnt++;
                totalSplitTris++;
            }
            else if(insertResult == VI_ENCROACHED)
            {
                /* If the newly inserted vertex encroaches upon a subsegment,
                 * delete the new vertex*/
                undoInsertVtxCDT();
            }
            else /* VI_VIOLATING*/
            {
                /* failed to insert the new vertex, but some subsegment was
                 * marked as being encroached */
            }
		}
		catch (...)
		{
            /* something unexpected happended*/
			failedSplitBadTris++;
            vertices.pop_back();
		}
	}

    VertexInsertType Triangulation::insertVtxCDT(IdxType v, const Triangle* tri)
	{
		PntTriLocationType locTri = locatePntTriangle(v, tri->V1(), tri->V2(), tri->V3());
		assert(locTri != PntTriLocationType::PT_OUTSIDE);
		if (locTri > PntTriLocationType::PT_INSIDE)  /*on edge*/
		{
			size_t edgeIdxInTri = static_cast<size_t>(locTri - PntTriLocationType::PT_ON_EDGE_1);
			size_t v1IdxInTri = edgeIdxInTri;
			size_t v2IdxInTri = (v1IdxInTri + 1) % 3;
			IdxType edgeV1 = tri->V(v1IdxInTri);
			IdxType edgeV2 = tri->V(v2IdxInTri);
			return insertVtxOnEdgeCDT(v, edgeV1, edgeV2);
		}
		else  /*inside triangle*/
		{
			return insertVtxInTriangleCDT(v, tri);
		}
	}

    void Triangulation::undoInsertVtxCDT()
    {
        
    }

    VertexInsertType Triangulation::insertVtxInTriangleCDT(IdxType v, const Triangle* tri)
	{
		IdxType v1 = tri->V1();
		IdxType v2 = tri->V2();
		IdxType v3 = tri->V3();
		deleteTriangle(v1, v2, v3);
		digCavityCDT(v, v1, v2);
		digCavityCDT(v, v2, v3);
		digCavityCDT(v, v3, v1);
	}

    VertexInsertType Triangulation::insertVtxOnEdgeCDT(IdxType v, IdxType edgeV1, IdxType edgeV2)
	{
		Edge edge(edgeV1, edgeV2);
		IdxType edgTriV = adjacent(edgeV1, edgeV2);
		IdxType edgRevTriV = adjacent(edgeV2, edgeV1);
		if (isVtxValid(edgTriV))
		{
			deleteTriangle(edgeV1, edgeV2, edgTriV);
			digCavityCDT(v, edgeV2, edgTriV);
			digCavityCDT(v, edgTriV, edgeV1);
		}
		if (isVtxValid(edgRevTriV))
		{
			deleteTriangle(edgeV2, edgeV1, edgRevTriV);
			digCavityCDT(v, edgeV1, edgRevTriV);
			digCavityCDT(v, edgRevTriV, edgeV2);
		}
        
        // mark two splitted edges constrained if it derived from a constrained edge
		if (isEdgeConstrained(edge) && (isVtxValid(edgTriV) || isVtxValid(edgRevTriV)))
		{
			removeConstrainEdge(edge);
			Edge splitEdge1(edgeV1, v);
			Edge splitEdge2(v, edgeV2);
			addConstrainEdge(splitEdge1);
			addConstrainEdge(splitEdge2);
		}
	}

	void Triangulation::digCavityCDT(IdxType v, IdxType v1, IdxType v2)
	{
		// do not flip constrained edge
		if (isEdgeConstrained(Edge(v1, v2)))
		{
			addTriangle(v, v1, v2);
			return;
		}

		IdxType vOpo = adjacent(v2, v1);
		if (!isVtxValid(vOpo))
		{
			addTriangle(v, v1, v2);
			return;
		}

		/* depth-first search to find all triangles that are no longer Delaunay*/
		if (isTriangleDelaunay(v, v1, v2, vOpo))
		{
			addTriangle(v, v1, v2);
		}
		else
		{
			deleteTriangle(v2, v1, vOpo);
			digCavityCDT(v, v1, vOpo);
			digCavityCDT(v, vOpo, v2);
		}
	}

	void Triangulation::locatePointWithGuide(IdxType v, const Triangle* guideTri, TriangulationLocation& loc)
	{
		/* Walking from a specific triangle to the triangle that contains v
		 *
		 * Walking is fast in practice if it follows two guides:
		 * 1.the vertices should be inserted in an order that has much spatial locality
		 * 2.each walk should begin at the most recently created triangle
		 *
		 * The typical walk visits small constant number of triangle
		 */

		assert(guideTri != nullptr);

		PntTriLocationType locTri = locatePntTriangle(v, guideTri->V1(), guideTri->V2(), guideTri->V3());
		if(locTri == PntTriLocationType::PT_INSIDE)
		{
			loc.type = TriangulationLocationType::TL_INSIDE;
			loc.tri = guideTri;
			return;
		}
		else if (locTri > PntTriLocationType::PT_INSIDE)  /*on edge*/
		{
			size_t edgeIdxInTri = static_cast<size_t>(locTri - PntTriLocationType::PT_ON_EDGE_1);
			size_t v1IdxInTri = edgeIdxInTri;
			size_t v2IdxInTri = (v1IdxInTri + 1) % 3;
			Edge edge(guideTri->V(v1IdxInTri), guideTri->V(v2IdxInTri));
			if (isEdgeConstrained(edge))
			{
				// do not insert vtx on constrained edge
				loc.type = TriangulationLocationType::TL_OUTSIDE;
				return;
			}
			loc.type = TriangulationLocationType::TL_INSIDE;
			loc.tri = guideTri;
			return;
		}

		PntLineLocationType plLocType;
		TriUSet visited;
		std::stack<const Triangle*> candidates;
		candidates.push(guideTri);
		const Triangle* neighborTri = nullptr;
		IdxType vStart;
		IdxType vEnd;
		while (!candidates.empty())
		{
			const Triangle* curTri = candidates.top();
			candidates.pop();
			assert(curTri != nullptr);

			bool foundInTriangle = true;
			// TODO: handle the case where vertex lie on the boundary
			for (IdxType i = 0; i < 3; ++i)
			{
				vStart = curTri->V(i % 3);
				vEnd = curTri->V((i + 1) % 3);
				plLocType = locatePntLine(v, vStart, vEnd);
				if (plLocType == PntLineLocationType::PL_RIGHT)
				{
					neighborTri = getTriangle(Edge(vEnd, vStart));
					if (neighborTri == nullptr)
					{
						foundInTriangle = false;
						break;
					}
					else if (visited.insert(*neighborTri).second)
					{
						foundInTriangle = false;
						candidates.push(neighborTri);
						break;
					}
				}
				else if (plLocType == PntLineLocationType::PL_ON_LINE)
				{
					Edge edge(vStart, vEnd);
					if (isEdgeConstrained(edge))
					{
						// do not insert vtx on constrained edge
						loc.type = TriangulationLocationType::TL_OUTSIDE;
						return;
					}
					foundInTriangle = true;
					break;
				}
			}

			if (foundInTriangle)
			{
				loc.type = TriangulationLocationType::TL_INSIDE;
				loc.tri = curTri;
				break;
			}
		}
	}

	bool Triangulation::isEdgeEncroached(const Edge& edge)
	{
		// check whether the vertex is in the diametral lens of the edge.
		IdxType v = adjacent(edge.V1(), edge.V2());
		if (isEdgeEncroachedByVtx(edge, v))
		{
			return true;
		}
		else
		{
			IdxType vOpo = adjacent(edge.V2(), edge.V1());
			if (vOpo == INVALID_VTX_IDX)
			{
				return false;
			}
			else
			{
				return isEdgeEncroachedByVtx(edge, vOpo);
			}
		}
	}

	bool Triangulation::isEdgeEncroachedByVtx(const Edge& edge, IdxType v)
	{
		/* A dot product of two sides of the triangle is used to check whether the
		 * angle at the vertex is greater than:
		 * 1. (90) degrees for diametral circles.
		 * 2. (180 - 2 `minangle') degrees for lenses.
		 */
		const Vertex& vtx1 = vertices[edge.V1()];
		const Vertex& vtx2 = vertices[edge.V2()];
		const Vertex& vtx = vertices[v];

		// (90) degrees for diametral circles
		PrecisionType dotproduct = (vtx1.X() - vtx.X()) * (vtx2.X() - vtx.X()) +
			(vtx1.Y() - vtx.Y()) * (vtx2.Y() - vtx.Y());

		if (dotproduct < 0)
		{
			// (180 - 2 `minangle') degrees for lenses
			PrecisionType a = dotproduct * dotproduct;
			PrecisionType b = (2.0 * goodAngleCosSquare - 1.0) * (2.0 * goodAngleCosSquare - 1.0);
			PrecisionType c = (vtx1.X() - vtx.X()) * (vtx1.X() - vtx.X()) +
				(vtx1.Y() - vtx.Y()) * (vtx1.Y() - vtx.Y());
			PrecisionType d = (vtx2.X() - vtx.X()) * (vtx2.X() - vtx.X()) +
				(vtx2.Y() - vtx.Y()) * (vtx2.Y() - vtx.Y());
			PrecisionType e = b * c * d;
			return a >= e;
		}
		return false;
	}

	void Triangulation::calTriangleQuality(const Triangle& tri, TriangleQuality& triQuality) const
	{
		Vertex tapex = vertices[tri.V1()];
		Vertex torg = vertices[tri.V2()];
		Vertex tdest = vertices[tri.V3()];
		PrecisionType dxod = torg.X() - tdest.X();
		PrecisionType dyod = torg.Y() - tdest.Y();
		PrecisionType dxda = tdest.X() - tapex.X();
		PrecisionType dyda = tdest.Y() - tapex.Y();
		PrecisionType dxao = tapex.X() - torg.X();
		PrecisionType dyao = tapex.Y() - torg.Y();
		PrecisionType dxod2 = dxod * dxod;
		PrecisionType dyod2 = dyod * dyod;
		PrecisionType dxda2 = dxda * dxda;
		PrecisionType dyda2 = dyda * dyda;
		PrecisionType dxao2 = dxao * dxao;
		PrecisionType dyao2 = dyao * dyao;

		// Find the lengths of the triangle's three edges
		PrecisionType apexlen = dxod2 + dyod2;
		PrecisionType orglen = dxda2 + dyda2;
		PrecisionType destlen = dxao2 + dyao2;

		PrecisionType minedge;
		PrecisionType minangle;
		if ((apexlen < orglen) && (apexlen < destlen))
		{
			// The edge opposite the apex is shortest
			minedge = apexlen;
			// Find the square of the cosine of the angle at the apex
			minangle = dxda * dxao + dyda * dyao;
			minangle = minangle * minangle / (orglen * destlen);
		}
		else if (orglen < destlen)
		{
			// The edge opposite the origin is shortest
			minedge = orglen;
			// Find the square of the cosine of the angle at the origin
			minangle = dxod * dxao + dyod * dyao;
			minangle = minangle * minangle / (apexlen * destlen);
		}
		else
		{
			// The edge opposite the destination is shortest
			minedge = destlen;
			// Find the square of the cosine of the angle at the destination
			minangle = dxod * dxda + dyod * dyda;
			minangle = minangle * minangle / (apexlen * orglen);
		}

		triQuality.area = 0.5 * (dxod * dyda - dyod * dxda);
		triQuality.minAngleCosSquare = minangle;
	}

	Vertex Triangulation::findCentriod(const Triangle& tri) const
	{
		const Vertex& v1 = vertices[tri.V1()];
		const Vertex& v2 = vertices[tri.V2()];
		const Vertex& v3 = vertices[tri.V3()];

		PrecisionType x = (v1.X() + v2.X() + v3.X()) / 3.0;
		PrecisionType y = (v1.Y() + v2.Y() + v3.Y()) / 3.0;

		return Vertex(x, y);
	}

	Vertex Triangulation::findCircumcenter(const Triangle& tri) const
	{
		const Vertex& v1 = vertices[tri.V1()];
		const Vertex& v2 = vertices[tri.V2()];
		const Vertex& v3 = vertices[tri.V3()];
		PrecisionType xdo = v3.X() - v2.X();
		PrecisionType ydo = v3.Y() - v2.Y();
		PrecisionType xao = v1.X() - v2.X();
		PrecisionType yao = v1.Y() - v2.Y();
		PrecisionType dodist = xdo * xdo + ydo * ydo;
		PrecisionType aodist = xao * xao + yao * yao;
		PrecisionType dadist = (v3.X() - v2.X()) * (v3.X() - v2.X()) +
			(v3.Y() - v2.Y()) * (v3.Y() - v2.Y());
		PrecisionType denominator = 0.5 / (xdo * yao - xao * ydo);
		PrecisionType dx = (yao * dodist - ydo * aodist) * denominator;
		PrecisionType dy = (xdo * aodist - xao * dodist) * denominator;

		/* Find the (squared) length of the triangle's shortest edge.  This
		 * serves as a conservative estimate of the insertion radius of the
		 * circumcenter's parent.  The estimate is used to ensure that
		 * the algorithm terminates even if very small angles appear in
		 * the input PSLG.*/
		if ((dodist < aodist) && (dodist < dadist))
		{
			// Find the position of the off-center, as described by Alper Ungor
			PrecisionType dxoff = 0.5 * xdo - offConstant * ydo;
			PrecisionType dyoff = 0.5 * ydo + offConstant * xdo;
			/* If the off-center is closer to the origin than the
			 * circumcenter, use the off-center instead*/
			if (dxoff * dxoff + dyoff * dyoff < dx * dx + dy * dy)
			{
				dx = dxoff;
				dy = dyoff;
			}
		}
		else if (aodist < dadist)
		{
			PrecisionType dxoff = 0.5 * xao + offConstant * yao;
			PrecisionType dyoff = 0.5 * yao - offConstant * xao;
			/* If the off-center is closer to the origin than the */
			/*   circumcenter, use the off-center instead.        */
			if (dxoff * dxoff + dyoff * dyoff < dx * dx + dy * dy)
			{
				dx = dxoff;
				dy = dyoff;
			}
		}
		else {
			PrecisionType dxoff = 0.5 * (v1.X() - v3.X()) -
				offConstant * (v1.Y() - v3.Y());
			PrecisionType dyoff = 0.5 * (v1.Y() - v3.Y()) +
				offConstant * (v1.X() - v3.X());
			/* If the off-center is closer to the destination than the */
			/*   circumcenter, use the off-center instead.             */
			if (dxoff * dxoff + dyoff * dyoff <
				(dx - xdo) * (dx - xdo) + (dy - ydo) * (dy - ydo)) {
				dx = xdo + dxoff;
				dy = ydo + dyoff;
			}
		}

		return Vertex(v2.X() + dx, v2.Y() + dy);
	}

	void Triangulation::addTriangle(IdxType v1, IdxType v2, IdxType v3)
	{
		assert(isTriangleValid(v1, v2, v3));

		Triangle tri(v1, v2, v3);
		triangles.push_back(tri);
		IdxType triIdx = triangles.size() - 1;

		// add each triangle for three times in 'edgeTriIdxTable'
		addEdgeTri(Edge(v1, v2), triIdx);
		addEdgeTri(Edge(v2, v3), triIdx);
		addEdgeTri(Edge(v3, v1), triIdx);
	}

	void Triangulation::deleteTriangle(IdxType v1, IdxType v2, IdxType v3)
	{
		assert(isTriangleValid(v1, v2, v3));

		// delete from 'edgeTriIdxTable'
		removeEdgeTri(Edge(v1, v2));
		removeEdgeTri(Edge(v2, v3));
		removeEdgeTri(Edge(v3, v1));
	}

	IdxType Triangulation::adjacent(IdxType v1, IdxType v2)
	{
		Edge edge(v1, v2);
		const Triangle* tri = getTriangle(edge);
		if (tri == nullptr)
			return INVALID_VTX_IDX;

		return tri->GetVtxCCW(v2);
	}

	std::pair<IdxType, IdxType> Triangulation::adjacentOne(IdxType v)
	{
		const Triangle* tri = nullptr;
		if (vtx2vtxMap.find(v) != vtx2vtxMap.end())
		{
			const IdxType& adjoinVtx = vtx2vtxMap[v];
			Edge edge(v, adjoinVtx);
			tri = getTriangle(edge);
		}

		if (tri == nullptr)
		{
			// try to search entire hash table
			for (auto& EdgeTriIter : edgeTriIdxTable)
			{
				const Edge& edge = EdgeTriIter.first;
				if (edge.IsContainVtx(v))
				{
					tri = getTriangle(edge);
					break;
				}
			}
		}

		assert(tri != nullptr);

		return std::make_pair(tri->GetVtxCCW(v), tri->GetVtxCW(v));
	}

	void Triangulation::addEdgeTri(const Edge& edge, IdxType triIdx)
	{
		if (isEdgeTriExist(edge))
		{
			// TODO: do some check here if edge already in 'edgeTriIdxTable'
			//assert(false);
		}

		edgeTriIdxTable.emplace(edge, triIdx);

		// cache first added triangle when inserting a vertex for the use of next vertex's point location
		if (mostRecentAddedTri == DUMMY_TRI_IDX)
		{
			mostRecentAddedTri = triIdx;
		}

		const Triangle* tri = getTriangle(triIdx);

		// cache added ghost triangles for the use of next vertex's point location
		if (isTriangleGhost(*tri))
		{
			Edge solidEdge(tri->GetVtxCCW(GHOST_VTX_IDX), tri->GetVtxCW(GHOST_VTX_IDX));
			ghostTriEdges.insert(solidEdge);
		}

		// cache most recent added triangle's vertices
		const IdxType& v1 = tri->V1();
		const IdxType& v2 = tri->V2();
		const IdxType& v3 = tri->V3();
		vtx2vtxMap[v1] = v2;
		vtx2vtxMap[v2] = v1;
		vtx2vtxMap[v3] = v2;
	}

	void Triangulation::removeEdgeTri(const Edge& edge)
	{
		if (!isEdgeTriExist(edge))
		{
			// TODO: do some check here if edge not in 'edgeTriIdxTable'
			//assert(false);
		}

		edgeTriIdxTable.erase(edge);
		ghostTriEdges.erase(edge);
	}

	bool Triangulation::isEdgeTriExist(const Edge& edge) const
	{
		return edgeTriIdxTable.find(edge) != edgeTriIdxTable.end();
	};

	bool Triangulation::isEdgeTriExist(const Triangle& tri) const
	{
		Edge triEdge1(tri.V1(), tri.V2());
		Edge triEdge2(tri.V2(), tri.V3());
		Edge triEdge3(tri.V3(), tri.V1());
		return isEdgeTriExist(triEdge1) &&
			isEdgeTriExist(triEdge2) &&
			isEdgeTriExist(triEdge3);
	};

	IdxType Triangulation::getTriangleIdx(const Edge& edge)
	{
		if (!isEdgeTriExist(edge)) return DUMMY_TRI_IDX;
		return edgeTriIdxTable[edge];
	};

	const Triangle* Triangulation::getTriangle(const Edge& edge)
	{
		IdxType triIdx = getTriangleIdx(edge);
		if (triIdx == DUMMY_TRI_IDX)
			return nullptr;
		return &triangles[triIdx];
	};

	const Triangle* Triangulation::getTriangle(IdxType triIdx)
	{
		if (triIdx == DUMMY_TRI_IDX)
			return nullptr;
		return &triangles[triIdx];
	};

	void Triangulation::addConstrainEdge(const Edge& edge)
	{
		constrainedEdges.insert(edge);
	}

	void Triangulation::removeConstrainEdge(const Edge& edge)
	{
		constrainedEdges.erase(edge);
		constrainedEdges.erase(Edge(edge.V2(), edge.V1()));
	}

	bool Triangulation::isEdgeConstrained(const Edge& edge) const
	{
		Edge reverseEdge(edge.V2(), edge.V1());
		return constrainedEdges.count(edge) || constrainedEdges.count(reverseEdge);
	}

	PntTriLocationType Triangulation::locatePntTriangle(IdxType v, IdxType v1, IdxType v2, IdxType v3) const
	{
		PntTriLocationType pntTriLoc = PntTriLocationType::PT_INSIDE;
		PntLineLocationType pntLineLoc = locatePntLine(v, v1, v2);
		if (pntLineLoc == PntLineLocationType::PL_RIGHT)
			return PntTriLocationType::PT_OUTSIDE;
		if (pntLineLoc == PntLineLocationType::PL_ON_LINE)
			pntTriLoc = PntTriLocationType::PT_ON_EDGE_1;
		pntLineLoc = locatePntLine(v, v2, v3);
		if (pntLineLoc == PntLineLocationType::PL_RIGHT)
			return PntTriLocationType::PT_OUTSIDE;
		if (pntLineLoc == PntLineLocationType::PL_ON_LINE)
			pntTriLoc = PntTriLocationType::PT_ON_EDGE_2;
		pntLineLoc = locatePntLine(v, v3, v1);
		if (pntLineLoc == PntLineLocationType::PL_RIGHT)
			return PntTriLocationType::PT_OUTSIDE;
		if (pntLineLoc == PntLineLocationType::PL_ON_LINE)
			pntTriLoc = PntTriLocationType::PT_ON_EDGE_3;
		return pntTriLoc;
	}

	PntLineLocationType Triangulation::locatePntLine(IdxType v, IdxType v1, IdxType v2) const
	{
		PredicateType orient = orient2dTest(v, v1, v2);
		if (orient < -orient2dTol)
			return PntLineLocationType::PL_RIGHT;
		if (orient > orient2dTol)
			return PntLineLocationType::PL_LEFT;
		return PntLineLocationType::PL_ON_LINE;
	}

	bool Triangulation::isPntInMidOfLine(IdxType v, IdxType v1, IdxType v2) const
	{
		/* NOTE: the vertex v must be collinear with line v1v2*/
		Vertex v1SubV = Vertex(
			vertices[v1].X() - vertices[v].X(),
			vertices[v1].Y() - vertices[v].Y()
		);
		Vertex v2SubV = Vertex(
			vertices[v2].X() - vertices[v].X(),
			vertices[v2].Y() - vertices[v].Y()
		);
		PrecisionType dotProduct = v1SubV.X() * v2SubV.X() + v1SubV.Y() * v2SubV.Y();
		return dotProduct < 0;
	}

	bool Triangulation::isVerticesCCW(IdxType v1, IdxType v2, IdxType v3) const
	{
		PredicateType orient2dTestValue = orient2dTest(v1, v2, v3);
		bool v1OnLeftOfV2V3 = orient2dTestValue > orient2dTol;
		return v1OnLeftOfV2V3;
	}

	bool Triangulation::isTriangleDelaunay(IdxType v1, IdxType v2, IdxType v3, IdxType v) const
	{
		PrecisionType incircleTestValue = incircleTest(v1, v2, v3, v);
		bool incircle = incircleTestValue > orient2dTol;
		return !incircle;
	}

	bool Triangulation::isTriangleValid(IdxType v1, IdxType v2, IdxType v3) const
	{
		if (!isVtxValid(v1) || !isVtxValid(v2) || !isVtxValid(v3))
			return false;

		if (!isTriangleGhost(v1, v2, v3) && !isVerticesCCW(v1, v2, v3))
			return false;

		return true;
	}

	PredicateType Triangulation::orient2dTest(IdxType v1, IdxType v2, IdxType v3) const
	{
		const Vertex& vtx1 = vertices[v1];
		const Vertex& vtx2 = vertices[v2];
		const Vertex& vtx3 = vertices[v3];
		return predicates::adaptive::orient2d(
			vtx1.X(), vtx1.Y(),
			vtx2.X(), vtx2.Y(),
			vtx3.X(), vtx3.Y()
		);
	}

	PredicateType Triangulation::incircleTest(IdxType v1, IdxType v2, IdxType v3, IdxType v) const
	{
		const Vertex& vtx1 = vertices[v1];
		const Vertex& vtx2 = vertices[v2];
		const Vertex& vtx3 = vertices[v3];
		const Vertex& vtx = vertices[v];
		PredicateType ret = predicates::adaptive::incircle(
			vtx1.X(), vtx1.Y(),
			vtx2.X(), vtx2.Y(),
			vtx3.X(), vtx3.Y(),
			vtx.X(), vtx.Y()
		);
		return ret;
	}

	IdxType Triangulation::getOppositeVtx(const Triangle& tri, IdxType v)
	{
		IdxType v1 = tri.GetVtxCCW(v);
		IdxType v2 = tri.GetVtxCW(v);
		Edge edge(v2, v1);
		const Triangle* oppositeTri = getTriangle(edge);
		if (oppositeTri == nullptr)
			return INVALID_VTX_IDX;
		return oppositeTri->GetVtxCCW(v1);
	}

	IdxType Triangulation::getOppositeTriangleIdx(const Triangle& tri, IdxType v)
	{
		Edge edge(tri.GetVtxCW(v), tri.GetVtxCCW(v));
		return getTriangleIdx(edge);
	}

	const Triangle* Triangulation::getOppositeTriangle(const Triangle& tri, IdxType v)
	{
		Edge edge(tri.GetVtxCW(v), tri.GetVtxCCW(v));
		return getTriangle(edge);
	}
}
