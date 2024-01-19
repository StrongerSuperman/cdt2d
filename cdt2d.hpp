#pragma once

#include <vector>
#include <stack>
#include <deque>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <limits>


namespace CDT
{
	#define CDT_ZERO 10E-9

	typedef size_t IdxType;
	typedef double PrecisionType;
	typedef unsigned short LayerDepth;
	typedef PrecisionType PredicateType;

	static IdxType INVALID_VTX_IDX = std::numeric_limits<IdxType>::max();
	static IdxType GHOST_VTX_IDX = std::numeric_limits<IdxType>::max() - 1;
	static IdxType DUMMY_TRI_IDX = std::numeric_limits<IdxType>::max();
	static LayerDepth DEFAULT_LAYER_DEPTH = std::numeric_limits<LayerDepth>::max();

	enum PntTriLocationType
	{
		PT_OUTSIDE = 0,
		PT_INSIDE,
		PT_ON_EDGE_1,
		PT_ON_EDGE_2,
		PT_ON_EDGE_3,
	};

	enum PntLineLocationType
	{
		PL_LEFT = 0,
		PL_RIGHT,
		PL_ON_LINE,
	};

	enum TriangulationLocationType
	{
		TL_UNKNOWN = 0,
		TL_OUTSIDE,
		TL_INSIDE,
	};

    enum VertexInsertType
    {
        VI_SUCCESS = 0,
        VI_ENCROACHED,
        VI_VIOLATING
    };

	class Vertex
	{
	public:
		explicit Vertex(PrecisionType x, PrecisionType y)
		{
			pos = std::make_pair(x, y);
		}
		inline PrecisionType X() const { return pos.first; };
		inline PrecisionType Y() const { return pos.second; };
	private:
		std::pair<PrecisionType, PrecisionType> pos;
        
	};

	class Edge
	{
	public:
		Edge()
		{
			vertices = std::make_pair(INVALID_VTX_IDX, INVALID_VTX_IDX);
		};
		Edge(IdxType v1, IdxType v2)
		{
			vertices = std::make_pair(v1, v2);
		};
		bool IsContainVtx(IdxType v) const
		{
			return vertices.first == v || vertices.second == v;
		};
		inline IdxType V1() const { return vertices.first; };
		inline IdxType V2() const { return vertices.second; };
		inline void SetV1(IdxType v) { vertices.first = v; };
		inline void SetV2(IdxType v) { vertices.second = v; };
		bool operator==(const Edge& edge) const
		{
			return vertices.first == edge.V1() && vertices.second == edge.V2();
		}
	private:
		std::pair<IdxType, IdxType> vertices;  // directed edge's vertices store in order
	};

	class Triangle
	{
	public:
		explicit Triangle(IdxType v1, IdxType v2, IdxType v3)
		{
			vertices.reserve(3);
			vertices.push_back(v1);
			vertices.push_back(v2);
			vertices.push_back(v3);
		};
		Triangle(const Triangle& tri)
		{
			vertices.resize(3);
			vertices[0] = tri.V1();
			vertices[1] = tri.V2();
			vertices[2] = tri.V3();
		}
		bool IsContainVtx(IdxType v) const
		{
			return vertices[0] == v || vertices[1] == v || vertices[2] == v;
		};
		IdxType GetVtxCCW(IdxType v) const
		{
			int i = 0;
			for (; i < 3; i++)
			{
				if (v == vertices[i])
					break;
			}
			return vertices[(i + 1) % 3];
		}
		IdxType GetVtxCW(IdxType v) const
		{
			int i = 0;
			for (; i < 3; i++)
			{
				if (v == vertices[i])
					break;
			}
			return vertices[(i + 2) % 3];
		}
		inline IdxType V1() const { return vertices[0]; };
		inline IdxType V2() const { return vertices[1]; };
		inline IdxType V3() const { return vertices[2]; };
		inline IdxType V(size_t idx) const { return vertices[idx]; };
		inline void SetV1(IdxType v) { vertices[0] = v; };
		inline void SetV2(IdxType v) { vertices[1] = v; };
		inline void SetV3(IdxType v) { vertices[2] = v; };
		bool operator==(const Triangle& tri) const
		{
			return (vertices[0] == tri.V1() && vertices[1] == tri.V2() && vertices[2] == tri.V3()) ||
				(vertices[1] == tri.V1() && vertices[2] == tri.V2() && vertices[0] == tri.V3()) || 
				(vertices[2] == tri.V1() && vertices[0] == tri.V2() && vertices[1] == tri.V3());
		}
	private:
		std::vector<IdxType> vertices;   // triangle's ccw vertices store in order
	};

	struct TriangleQuality
	{
		PrecisionType minAngleCosSquare;
		PrecisionType area;
	};

	class TriangulationLocation
	{
	public:
		TriangulationLocationType type;
		TriangulationLocation()
		{
			tri = nullptr;
			type = TriangulationLocationType::TL_UNKNOWN;
			onBoundaryEdge = false;
		}
		TriangulationLocation(const TriangulationLocation& rhs)
		{
			type = rhs.type;
			tri = rhs.tri;
			inHPEdges = rhs.inHPEdges;
			onHPEdges = rhs.onHPEdges;
			boundaryEdge = rhs.boundaryEdge;
			onBoundaryEdge = rhs.onBoundaryEdge;
		};
		~TriangulationLocation() {};
		const Triangle* tri;                // inside the triangle or on the triangle edge

		std::vector<Edge> inHPEdges;        // in the right half-plane of these directed edges
		std::vector<Edge> onHPEdges;        // on the right half-plane of these directed edges
		Edge boundaryEdge;                  // boundary edge
		bool onBoundaryEdge;                // whether on the boundary edge
	};

	inline static void hash_combine_value(std::size_t& seed, std::size_t hash_value)
	{
		seed ^= hash_value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
	}

	template <class T>
	inline static void hash_combine(std::size_t& seed, const T& v)
	{
		std::hash<T> hasher;
		hash_combine_value(seed, hasher(v));
	}

	struct HashEdge
	{
		size_t operator()(const Edge& edge) const
		{
			size_t value = 0;
			hash_combine(value, edge.V1());
			hash_combine(value, edge.V2());
			return value;
		}
	};

	struct HashTriangle
	{
		size_t operator()(const Triangle& tri) const
		{
			size_t value = 0;
			hash_combine(value, tri.V1());
			hash_combine(value, tri.V2());
			hash_combine(value, tri.V3());
			return value;
		}
	};

	typedef std::unordered_map<Edge, IdxType, HashEdge> EdgToTriIdxUMap;

	typedef std::unordered_set<Edge, HashEdge> EdgeUSet;
	typedef std::stack<Edge> EdgeStack;

	typedef std::unordered_set<Triangle, HashTriangle> TriUSet;
	typedef std::unordered_set<IdxType> TriIdxUSet;
	typedef std::stack<IdxType> TriIdxStack;
	typedef std::deque<IdxType> TriIdxDeque;

	typedef std::unordered_map<IdxType, LayerDepth> TriLayerUMap;
	typedef std::unordered_map<LayerDepth, TriIdxUSet> LayerTriIdxesUMap;

	class Triangulation
	{
	public:
		Triangulation();

		void InsertVertices(const std::vector<Vertex>& vertexList);
		void InsertEdges(const std::vector<Edge>& edgeList,
			const std::vector<std::vector<Edge>>& holesEdgeList = std::vector<std::vector<Edge>>());
		void Perform(bool delaunayRefine=true);

		const std::vector<Triangle>& GetTriangles() const { return triangles; };
		const std::vector<Vertex>& GetVertices() const { return vertices; };
		const EdgeUSet& GetConstrainedEdges() const { return constrainedEdges; };

	private:
		/* Insert vertex*/
		void insertVertex(IdxType v);
		void locatePoint(IdxType v, TriangulationLocation& loc);
		void insertVtxInsideTriangulation(IdxType v, const Triangle* tri);
		void insertVtxOnEdge(IdxType v, IdxType edgV1, IdxType edgV2);
		void insertVtxInTriangle(IdxType v, const Triangle* tri);
		void insertVtxOutsideTriangulation(IdxType v, const TriangulationLocation& loc);
		void digCavity(IdxType v, IdxType v1, IdxType v2);

		/* Insert edge*/
		void insertEdge(const Edge& edge);
		void insertEdgeIter(const Edge& edge, std::vector<Edge>& remainingEdges);
		void triangulatePseudopolygon(IdxType vStart, IdxType vEnd, std::vector<IdxType>&& midPts);
		size_t findDelaunayPoint(IdxType vStart, IdxType vEnd, const std::vector<IdxType>& midPts);

		/* Remove outers and holes*/
		void removeOuterAndHoles();
		void calculateTriangleDepths(std::vector<LayerDepth>& triDepths);
		TriLayerUMap peelLayer(TriIdxStack seeds, const LayerDepth layerDepth, std::vector<LayerDepth>& triDepths);

		/* Remap*/
		void remapTriangulation();

		/* Delaunay refine*/
		void delaunayRefinement();
		void splitEdge(const Edge& edge);
		void splitTriangle(IdxType badTriIdx);
        VertexInsertType insertVtxCDT(IdxType v, const Triangle* tri);
        void undoInsertVtxCDT();
        VertexInsertType insertVtxInTriangleCDT(IdxType v, const Triangle* tri);
        VertexInsertType insertVtxOnEdgeCDT(IdxType v, IdxType edgV1, IdxType edgV2);
		void digCavityCDT(IdxType v, IdxType v1, IdxType v2);
		void locatePointWithGuide(IdxType v, const Triangle* guideTri, TriangulationLocation& loc);
		bool isEdgeEncroached(const Edge& edge);
		bool isEdgeEncroachedByVtx(const Edge& edge, IdxType v);
		void calTriangleQuality(const Triangle& tri, TriangleQuality& triQuality) const;
		Vertex findCircumcenter(const Triangle& tri) const;
		Vertex findCentriod(const Triangle& tri) const;

	private:
		/* Add a positively oriented triangle v1v2v3, and do nothing if v1v2v3 is not
		positively oriented(ccw) or one edge-triangle exist in edgeTriIdxTable*/
		void addTriangle(IdxType v1, IdxType v2, IdxType v3);

		/* Delete a positively oriented(ccw) triangle v1v2v3*/
		void deleteTriangle(IdxType v1, IdxType v2, IdxType v3);

		/* Return a vertex v3 such that v1v2v3 is a positively oriented(ccw) triangle*/
		IdxType adjacent(IdxType v1, IdxType v2);

		/* Return a vertex v2, v3 such that v1v2v3 is a positively oriented(ccw) triangle*/
		std::pair<IdxType, IdxType> adjacentOne(IdxType v);

	private:
		void addEdgeTri(const Edge& edge, IdxType triIdx);
		void removeEdgeTri(const Edge& edge);
		bool isEdgeTriExist(const Edge& edge) const;
		bool isEdgeTriExist(const Triangle& tri) const;

		IdxType getTriangleIdx(const Edge& edge);
		const Triangle* getTriangle(const Edge& edge);
		const Triangle* getTriangle(IdxType triIdx);

		void addConstrainEdge(const Edge& edge);
		void removeConstrainEdge(const Edge& edge);
		bool isEdgeConstrained(const Edge& edge) const;

	private:
		PntTriLocationType locatePntTriangle(IdxType v, IdxType v1, IdxType v2, IdxType v3) const;
		PntLineLocationType locatePntLine(IdxType v, IdxType v1, IdxType v2) const;
		bool isPntInMidOfLine(IdxType v, IdxType v1, IdxType v2) const;

		PredicateType orient2dTest(IdxType v1, IdxType v2, IdxType v3) const;
		PredicateType incircleTest(IdxType v1, IdxType v2, IdxType v3, IdxType v) const;

		IdxType getOppositeVtx(const Triangle& tri, IdxType v);
		IdxType getOppositeTriangleIdx(const Triangle& tri, IdxType v);
		const Triangle* getOppositeTriangle(const Triangle& tri, IdxType v);

		bool isVerticesCCW(IdxType v1, IdxType v2, IdxType v3) const;
		bool isTriangleDelaunay(IdxType v1, IdxType v2, IdxType v3, IdxType v) const;
		bool isTriangleValid(IdxType v1, IdxType v2, IdxType v3) const;

		inline bool isTriangleGhost(const Triangle& tri) const
		{
			return isTriangleGhost(tri.V1(), tri.V2(), tri.V3());
		};
		inline bool isTriangleGhost(IdxType v1, IdxType v2, IdxType v3) const
		{
			return isVtxGhost(v1) || isVtxGhost(v2) || isVtxGhost(v3);
		};
		inline bool isEdgeGhost(const Edge& edge) const
		{
			return isEdgeGhost(edge.V1(), edge.V2());
		};
		inline bool isEdgeGhost(IdxType v1, IdxType v2) const
		{
			return isVtxGhost(v1) || isVtxGhost(v2);
		};
		inline bool isVtxGhost(IdxType v) const
		{
			return v == GHOST_VTX_IDX;
		};
		inline bool isVtxValid(IdxType v) const
		{
			return v != INVALID_VTX_IDX;
		};

	private:
		std::vector<Vertex> vertices;
		std::vector<Edge> edges;
		std::vector<std::vector<Edge>> holesEdges;
		std::vector<Triangle> triangles;

		/* each directed edge has a certain mapping triangle which can be solid/ghost triangle*/
		EdgToTriIdxUMap edgeTriIdxTable;
		/* TODO: use kd-tree to find nearest point*/
		IdxType mostRecentAddedTri = DUMMY_TRI_IDX;
		/* most recently added triangle adjoining u(key) also has a v(value) for a vertex*/
		std::unordered_map<IdxType, IdxType> vtx2vtxMap;    
		/* a unordered set to store all ghost edges*/
		EdgeUSet ghostTriEdges;
		/* a unordered set to store all constrained edges*/
		EdgeUSet constrainedEdges;
        
        /* for delaunay refinement*/
        EdgeStack encroachedEdges;
        TriIdxDeque badTriangles;

		PredicateType orient2dTol = 0;
		PredicateType incircleTol = 0;

		int insertSteinerCnt = 0;              // current insert steiner point cnt
		int maxSteinerCnt = 500;               // max insert steiner pnt cnt
		PrecisionType maxArea = 10.0;          // maximum area bound
		PrecisionType minAngle = 25.0;         // minimum angle bound

		PrecisionType goodAngleCosSquare;      // cosine squared of minAngle
		PrecisionType offConstant;             // constant used to place off-center Steiner points

		/* for debug*/
		int failedSplitBadTris = 0;
		int totalSplitTris = 0;
		int totalSplitEdges = 0;
	};
}
