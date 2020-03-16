#include"RTree.H"

#define RTREE_TEMPLATE template<class DATATYPE, class ELEMTYPE, int NUMDIMS, class ELEMTYPEREAL, int TMAXNODES, int TMINNODES>
#define RTREE_QUAL RTree<DATATYPE  , ELEMTYPE, NUMDIMS, ELEMTYPEREAL, TMAXNODES, TMINNODES>

RTREE_TEMPLATE
RTREE_QUAL::RTree()
{
	ASSERT(MAXNODES > MINNODES);
	ASSERT(MINNODES > 0);


	// We only support machine word size simple data type eg. integer index or object pointer.
	// Since we are storing as union with non data branch
	//ASSERT(sizeof(DATATYPE) == sizeof(void*) || sizeof(DATATYPE) == sizeof(int));

	// Precomputed volumes of the unit spheres for the first few dimensions
	const float UNIT_SPHERE_VOLUMES[] =
	{
		0.000000f, 2.000000f, 3.141593f, // Dimension  0,1,2
		4.188790f, 4.934802f, 5.263789f, // Dimension  3,4,5
		5.167713f, 4.724766f, 4.058712f, // Dimension  6,7,8
		3.298509f, 2.550164f, 1.884104f, // Dimension  9,10,11
		1.335263f, 0.910629f, 0.599265f, // Dimension  12,13,14
		0.381443f, 0.235331f, 0.140981f, // Dimension  15,16,17
		0.082146f, 0.046622f, 0.025807f, // Dimension  18,19,20
	};

	m_root = AllocNode();
	m_root->m_level = 0;
	m_unitSphereVolume =static_cast<ELEMTYPEREAL>(UNIT_SPHERE_VOLUMES[NUMDIMS]);
}


RTREE_TEMPLATE
RTREE_QUAL::~RTree()
{
	Reset(); // Free, or reset node memory
}


RTREE_TEMPLATE
void RTREE_QUAL::Insert(const ELEMTYPE a_min[NUMDIMS], const ELEMTYPE a_max[NUMDIMS], const DATATYPE& a_dataId)
{
#ifdef _DEBUG
	for(int index=0; index<NUMDIMS; ++index)
	{
		ASSERT(a_min[index] <= a_max[index]);
	}
#endif //_DEBUG

	Rect rect;

	for(int axis=0; axis<NUMDIMS; ++axis)
	{
		rect.m_min[axis] = a_min[axis];
		rect.m_max[axis] = a_max[axis];
	}

	InsertRect(&rect, a_dataId, &m_root, 0);
}


RTREE_TEMPLATE
void RTREE_QUAL::Remove(const ELEMTYPE a_min[NUMDIMS], const ELEMTYPE a_max[NUMDIMS], const DATATYPE& a_dataId)
{
#ifdef _DEBUG
	for(int index=0; index<NUMDIMS; ++index)
	{
		ASSERT(a_min[index] <= a_max[index]);
	}
#endif //_DEBUG

	Rect rect;

	for(int axis=0; axis<NUMDIMS; ++axis)
	{
		rect.m_min[axis] = a_min[axis];
		rect.m_max[axis] = a_max[axis];
	}

	RemoveRect(&rect, a_dataId, &m_root);
}

RTREE_TEMPLATE
int RTREE_QUAL::Search(const ELEMTYPE a_x, const ELEMTYPE a_y, bool/* __cdecl*/ a_resultCallback(DATATYPE a_data, void* a_context), void* a_context)
{
#ifdef _DEBUG
	for(int index=0; index<NUMDIMS; ++index)
	{
		ASSERT(a_min[index] <= a_max[index]);
	}
#endif //_DEBUG

	Rect rect;


	rect.m_min[0] = a_x;
	rect.m_max[0] = a_x;
	rect.m_min[1] = a_y;
	rect.m_max[1] = a_y;


	// NOTE: May want to return search result another way, perhaps returning the number of found elements here.

	int foundCount = 0;
	Search(m_root, &rect, foundCount, a_context, a_resultCallback);

	return foundCount;
}



RTREE_TEMPLATE
int RTREE_QUAL::Search(const ELEMTYPE a_min[NUMDIMS], const ELEMTYPE a_max[NUMDIMS], bool/* __cdecl*/  a_resultCallback(DATATYPE a_data, void* a_context), void* a_context)
{
#ifdef _DEBUG
	for(int index=0; index<NUMDIMS; ++index)
	{
		ASSERT(a_min[index] <= a_max[index]);
	}
#endif //_DEBUG

	Rect rect;

	for(int axis=0; axis<NUMDIMS; ++axis)
	{
		rect.m_min[axis] = a_min[axis];
		rect.m_max[axis] = a_max[axis];
	}

	// NOTE: May want to return search result another way, perhaps returning the number of found elements here.

	int foundCount = 0;
	Search(m_root, &rect, foundCount, a_resultCallback, a_context);

	return foundCount;
}

RTREE_TEMPLATE
int RTREE_QUAL::Search(const ELEMTYPE a[NUMDIMS], void* a_context, bool/* __cdecl*/  a_resultCallback(DATATYPE a_data, void* a_context))
{

	Point point;

	for(int axis=0; axis<NUMDIMS; ++axis)
		point.m[axis] = a[axis];

	// NOTE: May want to return search result another way, perhaps returning the number of found elements here.
 
	int foundCount = 0;
	Search(m_root, &point, foundCount, a_context, a_resultCallback);

	return foundCount;
}
RTREE_TEMPLATE
bool RTREE_QUAL::returndata(DATATYPE a_data, void* a_context)
{
	DATATYPE * pdata=(reinterpret_cast<DATATYPE *> (a_context));
	*pdata=a_data;


	return false;//stop
	return true; // keep going
}


RTREE_TEMPLATE
int RTREE_QUAL::Count()
{
	int count = 0;
	CountRec(m_root, count);

	return count;
}



RTREE_TEMPLATE
void RTREE_QUAL::CountRec(Node* a_node, int& a_count)
{
	if(a_node->IsInternalNode())   // not a leaf node
	{
		for(int index = 0; index < a_node->m_count; ++index)
		{
			CountRec(a_node->m_branch[index].m_child, a_count);
		}
	}
	else     // A leaf node
	{
		a_count += a_node->m_count;
	}
}


RTREE_TEMPLATE
bool RTREE_QUAL::Load(const char* a_fileName)
{
	RemoveAll(); // Clear existing tree

	RTFileStream stream;
	if(!stream.OpenRead(a_fileName))
	{
		return false;
	}

	bool result = Load(stream);

	stream.Close();

	return result;
};



RTREE_TEMPLATE
bool RTREE_QUAL::Load(RTFileStream& a_stream)
{
	// Write some kind of header
	int _dataFileId = ('R'<<0)|('T'<<8)|('R'<<16)|('E'<<24);
	int _dataSize = sizeof(DATATYPE);
	int _dataNumDims = NUMDIMS;
	int _dataElemSize = sizeof(ELEMTYPE);
	int _dataElemRealSize = sizeof(ELEMTYPEREAL);
	int _dataMaxNodes = TMAXNODES;
	int _dataMinNodes = TMINNODES;

	int dataFileId = 0;
	int dataSize = 0;
	int dataNumDims = 0;
	int dataElemSize = 0;
	int dataElemRealSize = 0;
	int dataMaxNodes = 0;
	int dataMinNodes = 0;

	a_stream.Read(dataFileId);
	a_stream.Read(dataSize);
	a_stream.Read(dataNumDims);
	a_stream.Read(dataElemSize);
	a_stream.Read(dataElemRealSize);
	a_stream.Read(dataMaxNodes);
	a_stream.Read(dataMinNodes);

	bool result = false;

	// Test if header was valid and compatible
	if(    (dataFileId == _dataFileId)
	        && (dataSize == _dataSize)
	        && (dataNumDims == _dataNumDims)
	        && (dataElemSize == _dataElemSize)
	        && (dataElemRealSize == _dataElemRealSize)
	        && (dataMaxNodes == _dataMaxNodes)
	        && (dataMinNodes == _dataMinNodes)
	  )
	{
		// Recursively load tree
		result = LoadRec(m_root, a_stream);
	}

	return result;
}


RTREE_TEMPLATE
bool RTREE_QUAL::LoadRec(Node* a_node, RTFileStream& a_stream)
{
	a_stream.Read(a_node->m_level);
	a_stream.Read(a_node->m_count);

	if(a_node->IsInternalNode())   // not a leaf node
	{
		for(int index = 0; index < a_node->m_count; ++index)
		{
			Branch* curBranch = &a_node->m_branch[index];

			a_stream.ReadArray(curBranch->m_rect.m_min, NUMDIMS);
			a_stream.ReadArray(curBranch->m_rect.m_max, NUMDIMS);

			curBranch->m_child = AllocNode();
			LoadRec(curBranch->m_child, a_stream);
		}
	}
	else     // A leaf node
	{
		for(int index = 0; index < a_node->m_count; ++index)
		{
			Branch* curBranch = &a_node->m_branch[index];

			a_stream.ReadArray(curBranch->m_rect.m_min, NUMDIMS);
			a_stream.ReadArray(curBranch->m_rect.m_max, NUMDIMS);

			a_stream.Read(curBranch->m_data);
		}
	}

	return true; // Should do more error checking on I/O operations
}


RTREE_TEMPLATE
bool RTREE_QUAL::Save(const char* a_fileName)
{
	RTFileStream stream;
	if(!stream.OpenWrite(a_fileName))
	{
		return false;
	}

	bool result = Save(stream);

	stream.Close();

	return result;
}


RTREE_TEMPLATE
bool RTREE_QUAL::Save(RTFileStream& a_stream)
{
	// Write some kind of header
	int dataFileId = ('R'<<0)|('T'<<8)|('R'<<16)|('E'<<24);
	int dataSize = sizeof(DATATYPE);
	int dataNumDims = NUMDIMS;
	int dataElemSize = sizeof(ELEMTYPE);
	int dataElemRealSize = sizeof(ELEMTYPEREAL);
	int dataMaxNodes = TMAXNODES;
	int dataMinNodes = TMINNODES;

	a_stream.Write(dataFileId);
	a_stream.Write(dataSize);
	a_stream.Write(dataNumDims);
	a_stream.Write(dataElemSize);
	a_stream.Write(dataElemRealSize);
	a_stream.Write(dataMaxNodes);
	a_stream.Write(dataMinNodes);

	// Recursively save tree
	bool result = SaveRec(m_root, a_stream);

	return result;
}


RTREE_TEMPLATE
bool RTREE_QUAL::SaveRec(Node* a_node, RTFileStream& a_stream)
{
	a_stream.Write(a_node->m_level);
	a_stream.Write(a_node->m_count);

	if(a_node->IsInternalNode())   // not a leaf node
	{
		for(int index = 0; index < a_node->m_count; ++index)
		{
			Branch* curBranch = &a_node->m_branch[index];

			a_stream.WriteArray(curBranch->m_rect.m_min, NUMDIMS);
			a_stream.WriteArray(curBranch->m_rect.m_max, NUMDIMS);

			SaveRec(curBranch->m_child, a_stream);
		}
	}
	else     // A leaf node
	{
		for(int index = 0; index < a_node->m_count; ++index)
		{
			Branch* curBranch = &a_node->m_branch[index];

			a_stream.WriteArray(curBranch->m_rect.m_min, NUMDIMS);
			a_stream.WriteArray(curBranch->m_rect.m_max, NUMDIMS);

			a_stream.Write(curBranch->m_data);
		}
	}

	return true; // Should do more error checking on I/O operations
}


RTREE_TEMPLATE
void RTREE_QUAL::RemoveAll()
{
	// Delete all existing nodes
	Reset();

	m_root = AllocNode();
	m_root->m_level = 0;
}


RTREE_TEMPLATE
void RTREE_QUAL::Reset()
{
#ifdef RTREE_DONT_USE_MEMPOOLS
	// Delete all existing nodes
	RemoveAllRec(m_root);
#else // RTREE_DONT_USE_MEMPOOLS
	// Just reset memory pools.  We are not using complex types
	// EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


RTREE_TEMPLATE
void RTREE_QUAL::RemoveAllRec(Node* a_node)
{
	ASSERT(a_node);
	ASSERT(a_node->m_level >= 0);

	if(a_node->IsInternalNode())   // This is an internal node in the tree
	{
		for(int index=0; index < a_node->m_count; ++index)
		{
			RemoveAllRec(a_node->m_branch[index].m_child);
		}
	}
	FreeNode(a_node);
}


RTREE_TEMPLATE
typename RTREE_QUAL::Node* RTREE_QUAL::AllocNode()
{
	Node* newNode;
#ifdef RTREE_DONT_USE_MEMPOOLS
	newNode = new Node;
#else // RTREE_DONT_USE_MEMPOOLS
	// EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
	InitNode(newNode);
	return newNode;
}


RTREE_TEMPLATE
void RTREE_QUAL::FreeNode(Node* a_node)
{
	ASSERT(a_node);

#ifdef RTREE_DONT_USE_MEMPOOLS
	delete a_node;
#else // RTREE_DONT_USE_MEMPOOLS
	// EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


// Allocate space for a node in the list used in DeletRect to
// store Nodes that are too empty.
RTREE_TEMPLATE
typename RTREE_QUAL::ListNode* RTREE_QUAL::AllocListNode()
{
#ifdef RTREE_DONT_USE_MEMPOOLS
	return new ListNode;
#else // RTREE_DONT_USE_MEMPOOLS
	// EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


RTREE_TEMPLATE
void RTREE_QUAL::FreeListNode(ListNode* a_listNode)
{
#ifdef RTREE_DONT_USE_MEMPOOLS
	delete a_listNode;
#else // RTREE_DONT_USE_MEMPOOLS
	// EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


RTREE_TEMPLATE
void RTREE_QUAL::InitNode(Node* a_node)
{
	a_node->m_count = 0;
	a_node->m_level = -1;
}


RTREE_TEMPLATE
void RTREE_QUAL::InitRect(Rect* a_rect)
{
	for(int index = 0; index < NUMDIMS; ++index)
	{
		a_rect->m_min[index] = static_cast<ELEMTYPE>(0);
		a_rect->m_max[index] = static_cast<ELEMTYPE>(0);
	}
}


// Inserts a new data rectangle into the index structure.
// Recursively descends tree, propagates splits back up.
// Returns 0 if node was not split.  Old node updated.
// If node was split, returns 1 and sets the pointer pointed to by
// new_node to point to the new node.  Old node updated to become one of two.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
RTREE_TEMPLATE
bool RTREE_QUAL::InsertRectRec(Rect* a_rect, const DATATYPE& a_id, Node* a_node, Node** a_newNode, int a_level)
{
	ASSERT(a_rect && a_node && a_newNode);
	ASSERT(a_level >= 0 && a_level <= a_node->m_level);

	int index;
	Branch branch;
	Node* otherNode;

	// Still above level for insertion, go down tree recursively
	if(a_node->m_level > a_level)
	{
		index = PickBranch(a_rect, a_node);
		if (!InsertRectRec(a_rect, a_id, a_node->m_branch[index].m_child, &otherNode, a_level))
		{
			// Child was not split
			a_node->m_branch[index].m_rect = CombineRect(a_rect, &(a_node->m_branch[index].m_rect));
			return false;
		}
		else     // Child was split
		{
			a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
			branch.m_child = otherNode;
			branch.m_rect = NodeCover(otherNode);
			return AddBranch(&branch, a_node, a_newNode);
		}
	}
	else if(a_node->m_level == a_level)     // Have reached level for insertion. Add rect, split if necessary
	{
		branch.m_rect = *a_rect;
		//branch.m_child = (Node*) a_id;
		branch.m_data =  a_id;
		// Child field of leaves contains id of data record
		return AddBranch(&branch, a_node, a_newNode);
	}
	else
	{
		// Should never occur
		ASSERT(0);
		return false;
	}
}


// Insert a data rectangle into an index structure.
// InsertRect provides for splitting the root;
// returns 1 if root was split, 0 if it was not.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
// InsertRect2 does the recursion.
//
RTREE_TEMPLATE
bool RTREE_QUAL::InsertRect(Rect* a_rect, const DATATYPE& a_id, Node** a_root, int a_level)
{
	ASSERT(a_rect && a_root);
	ASSERT(a_level >= 0 && a_level <= (*a_root)->m_level);
#ifdef _DEBUG
	for(int index=0; index < NUMDIMS; ++index)
	{
		ASSERT(a_rect->m_min[index] <= a_rect->m_max[index]);
	}
#endif //_DEBUG  

	Node* newRoot;
	Node* newNode;
	Branch branch;

	if(InsertRectRec(a_rect, a_id, *a_root, &newNode, a_level))   // Root split
	{
		newRoot = AllocNode();  // Grow tree taller and new root
		newRoot->m_level = (*a_root)->m_level + 1;
		branch.m_rect = NodeCover(*a_root);
		branch.m_child = *a_root;
		AddBranch(&branch, newRoot, NULL);
		branch.m_rect = NodeCover(newNode);
		branch.m_child = newNode;
		AddBranch(&branch, newRoot, NULL);
		*a_root = newRoot;
		return true;
	}

	return false;
}


// Find the smallest rectangle that includes all rectangles in branches of a node.
RTREE_TEMPLATE
typename RTREE_QUAL::Rect RTREE_QUAL::NodeCover(Node* a_node)
{
	ASSERT(a_node);

	int firstTime = true;
	Rect rect;
	InitRect(&rect);

	for(int index = 0; index < a_node->m_count; ++index)
	{
		if(firstTime)
		{
			rect = a_node->m_branch[index].m_rect;
			firstTime = false;
		}
		else
		{
			rect = CombineRect(&rect, &(a_node->m_branch[index].m_rect));
		}
	}

	return rect;
}


// Add a branch to a node.  Split the node if necessary.
// Returns 0 if node not split.  Old node updated.
// Returns 1 if node split, sets *new_node to address of new node.
// Old node updated, becomes one of two.
RTREE_TEMPLATE
bool RTREE_QUAL::AddBranch(Branch* a_branch, Node* a_node, Node** a_newNode)
{
	ASSERT(a_branch);
	ASSERT(a_node);

	if(a_node->m_count < MAXNODES)   // Split won't be necessary
	{
		a_node->m_branch[a_node->m_count] = *a_branch;
		++a_node->m_count;

		return false;
	}
	else
	{
		ASSERT(a_newNode);

		SplitNode(a_node, a_branch, a_newNode);
		return true;
	}
}


// Disconnect a dependent node.
// Caller must return (or stop using iteration index) after this as count has changed
RTREE_TEMPLATE
void RTREE_QUAL::DisconnectBranch(Node* a_node, int a_index)
{
	ASSERT(a_node && (a_index >= 0) && (a_index < MAXNODES));
	ASSERT(a_node->m_count > 0);

	// Remove element by swapping with the last element to prevent gaps in array
	a_node->m_branch[a_index] = a_node->m_branch[a_node->m_count - 1];

	--a_node->m_count;
}


// Pick a branch.  Pick the one that will need the smallest increase
// in area to accomodate the new rectangle.  This will result in the
// least total area for the covering rectangles in the current node.
// In case of a tie, pick the one which was smaller before, to get
// the best resolution when searching.
RTREE_TEMPLATE
int RTREE_QUAL::PickBranch(Rect* a_rect, Node* a_node)
{
	ASSERT(a_rect && a_node);

	bool firstTime = true;
	ELEMTYPEREAL increase;
	ELEMTYPEREAL bestIncr =static_cast<ELEMTYPEREAL>(-1);
	ELEMTYPEREAL area;
	ELEMTYPEREAL bestArea;
	int best=-1;
	Rect tempRect;

	for(int index=0; index < a_node->m_count; ++index)
	{
		Rect* curRect = &a_node->m_branch[index].m_rect;
		area = CalcRectVolume(curRect);
		tempRect = CombineRect(a_rect, curRect);
		increase = CalcRectVolume(&tempRect) - area;
		if((increase < bestIncr) || firstTime)
		{
			best = index;
			bestArea = area;
			bestIncr = increase;
			firstTime = false;
		}
		else if((increase == bestIncr) && (area < bestArea))
		{
			best = index;
			bestArea = area;
			bestIncr = increase;
		}
	}
	return best;
}


// Combine two rectangles into larger one containing both
RTREE_TEMPLATE
typename RTREE_QUAL::Rect RTREE_QUAL::CombineRect(Rect* a_rectA, Rect* a_rectB)
{
	ASSERT(a_rectA && a_rectB);

	Rect newRect;

	for(int index = 0; index < NUMDIMS; ++index)
	{
		newRect.m_min[index] = Min(a_rectA->m_min[index], a_rectB->m_min[index]);
		newRect.m_max[index] = Max(a_rectA->m_max[index], a_rectB->m_max[index]);
	}

	return newRect;
}



// Split a node.
// Divides the nodes branches and the extra one between two nodes.
// Old node is one of the new ones, and one really new one is created.
// Tries more than one method for choosing a partition, uses best result.
RTREE_TEMPLATE
void RTREE_QUAL::SplitNode(Node* a_node, Branch* a_branch, Node** a_newNode)
{
	ASSERT(a_node);
	ASSERT(a_branch);

	// Could just use local here, but member or external is faster since it is reused
	PartitionVars localVars;
	PartitionVars* parVars = &localVars;
	int level;

	// Load all the branches into a buffer, initialize old node
	level = a_node->m_level;
	GetBranches(a_node, a_branch, parVars);

	// Find partition
	ChoosePartition(parVars, MINNODES);

	// Put branches from buffer into 2 nodes according to chosen partition
	*a_newNode = AllocNode();
	(*a_newNode)->m_level = a_node->m_level = level;
	LoadNodes(a_node, *a_newNode, parVars);

	ASSERT((a_node->m_count + (*a_newNode)->m_count) == parVars->m_total);
}


// Calculate the n-dimensional volume of a rectangle
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_QUAL::RectVolume(Rect* a_rect)
{
	ASSERT(a_rect);

	ELEMTYPEREAL volume = static_cast<ELEMTYPEREAL>(1);

	for(int index=0; index<NUMDIMS; ++index)
	{
		volume *= a_rect->m_max[index] - a_rect->m_min[index];
	}

	ASSERT(volume >= static_cast<ELEMTYPEREAL>(0));

	return volume;
}


// The exact volume of the bounding sphere for the given Rect
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_QUAL::RectSphericalVolume(Rect* a_rect)
{
	ASSERT(a_rect);

	ELEMTYPEREAL sumOfSquares =static_cast<ELEMTYPEREAL>(0) ;
	ELEMTYPEREAL radius;

	for(int index=0; index < NUMDIMS; ++index)
	{
		ELEMTYPEREAL halfExtent = ((static_cast<ELEMTYPEREAL>(a_rect->m_max[index]) - static_cast<ELEMTYPEREAL>(a_rect->m_min[index]) )* 0.5f);
		sumOfSquares += halfExtent * halfExtent;
	}

	radius =static_cast<ELEMTYPEREAL>(sqrt(sumOfSquares));

	// Pow maybe slow, so test for common dims like 2,3 and just use x*x, x*x*x.
	if(NUMDIMS == 3)
	{
		return (radius * radius * radius * m_unitSphereVolume);
	}
	else if(NUMDIMS == 2)
	{
		return (radius * radius * m_unitSphereVolume);
	}
	else
	{
		return static_cast<ELEMTYPEREAL>(pow(radius, NUMDIMS) * m_unitSphereVolume);
	}
}


// Use one of the methods to calculate retangle volume
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_QUAL::CalcRectVolume(Rect* a_rect)
{
#ifdef RTREE_USE_SPHERICAL_VOLUME
	return RectSphericalVolume(a_rect); // Slower but helps certain merge cases
#else // RTREE_USE_SPHERICAL_VOLUME
	return RectVolume(a_rect); // Faster but can cause poor merges
#endif // RTREE_USE_SPHERICAL_VOLUME  
}


// Load branch buffer with branches from full node plus the extra branch.
RTREE_TEMPLATE
void RTREE_QUAL::GetBranches(Node* a_node, Branch* a_branch, PartitionVars* a_parVars)
{
	ASSERT(a_node);
	ASSERT(a_branch);

	ASSERT(a_node->m_count == MAXNODES);

	// Load the branch buffer
	for(int index=0; index < MAXNODES; ++index)
	{
		a_parVars->m_branchBuf[index] = a_node->m_branch[index];
	}
	a_parVars->m_branchBuf[MAXNODES] = *a_branch;
	a_parVars->m_branchCount = MAXNODES + 1;

	// Calculate rect containing all in the set
	a_parVars->m_coverSplit = a_parVars->m_branchBuf[0].m_rect;
	for(int index=1; index < MAXNODES+1; ++index)
	{
		a_parVars->m_coverSplit = CombineRect(&a_parVars->m_coverSplit, &a_parVars->m_branchBuf[index].m_rect);
	}
	a_parVars->m_coverSplitArea = CalcRectVolume(&a_parVars->m_coverSplit);

	InitNode(a_node);
}


// Method #0 for choosing a partition:
// As the seeds for the two groups, pick the two rects that would waste the
// most area if covered by a single rectangle, i.e. evidently the worst pair
// to have in the same group.
// Of the remaining, one at a time is chosen to be put in one of the two groups.
// The one chosen is the one with the greatest difference in area expansion
// depending on which group - the rect most strongly attracted to one group
// and repelled from the other.
// If one group gets too full (more would force other group to violate min
// fill requirement) then other group gets the rest.
// These last are the ones that can go in either group most easily.
RTREE_TEMPLATE
void RTREE_QUAL::ChoosePartition(PartitionVars* a_parVars, int a_minFill)
{
	ASSERT(a_parVars);

	ELEMTYPEREAL biggestDiff;
	int group, chosen=-1, betterGroup=-1;

	InitParVars(a_parVars, a_parVars->m_branchCount, a_minFill);
	PickSeeds(a_parVars);

	while (((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total)
	        && (a_parVars->m_count[0] < (a_parVars->m_total - a_parVars->m_minFill))
	        && (a_parVars->m_count[1] < (a_parVars->m_total - a_parVars->m_minFill)))
	{
		biggestDiff = static_cast<ELEMTYPEREAL>(-1);
		for(int index=0; index<a_parVars->m_total; ++index)
		{
			if(!a_parVars->m_taken[index])
			{
				Rect* curRect = &a_parVars->m_branchBuf[index].m_rect;
				Rect rect0 = CombineRect(curRect, &a_parVars->m_cover[0]);
				Rect rect1 = CombineRect(curRect, &a_parVars->m_cover[1]);
				ELEMTYPEREAL growth0 = CalcRectVolume(&rect0) - a_parVars->m_area[0];
				ELEMTYPEREAL growth1 = CalcRectVolume(&rect1) - a_parVars->m_area[1];
				ELEMTYPEREAL diff = growth1 - growth0;
				if(diff >= 0)
				{
					group = 0;
				}
				else
				{
					group = 1;
					diff = -diff;
				}

				if(diff > biggestDiff)
				{
					biggestDiff = diff;
					chosen = index;
					betterGroup = group;
				}
				else if((diff == biggestDiff) && (a_parVars->m_count[group] < a_parVars->m_count[betterGroup]))
				{
					chosen = index;
					betterGroup = group;
				}
			}
		}
		Classify(chosen, betterGroup, a_parVars);
	}

	// If one group too full, put remaining rects in the other
	if((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total)
	{
		if(a_parVars->m_count[0] >= a_parVars->m_total - a_parVars->m_minFill)
		{
			group = 1;
		}
		else
		{
			group = 0;
		}
		for(int index=0; index<a_parVars->m_total; ++index)
		{
			if(!a_parVars->m_taken[index])
			{
				Classify(index, group, a_parVars);
			}
		}
	}

	ASSERT((a_parVars->m_count[0] + a_parVars->m_count[1]) == a_parVars->m_total);
	ASSERT((a_parVars->m_count[0] >= a_parVars->m_minFill) &&
	       (a_parVars->m_count[1] >= a_parVars->m_minFill));
}


// Copy branches from the buffer into two nodes according to the partition.
RTREE_TEMPLATE
void RTREE_QUAL::LoadNodes(Node* a_nodeA, Node* a_nodeB, PartitionVars* a_parVars)
{
	ASSERT(a_nodeA);
	ASSERT(a_nodeB);
	ASSERT(a_parVars);

	for(int index=0; index < a_parVars->m_total; ++index)
	{
		ASSERT(a_parVars->m_partition[index] == 0 || a_parVars->m_partition[index] == 1);

		if(a_parVars->m_partition[index] == 0)
		{
			AddBranch(&a_parVars->m_branchBuf[index], a_nodeA, NULL);
		}
		else if(a_parVars->m_partition[index] == 1)
		{
			AddBranch(&a_parVars->m_branchBuf[index], a_nodeB, NULL);
		}
	}
}


// Initialize a PartitionVars structure.
RTREE_TEMPLATE
void RTREE_QUAL::InitParVars(PartitionVars* a_parVars, int a_maxRects, int a_minFill)
{
	ASSERT(a_parVars);

	a_parVars->m_count[0] = a_parVars->m_count[1] = 0;
	a_parVars->m_area[0] = a_parVars->m_area[1] =static_cast<ELEMTYPEREAL>(0) ;
	a_parVars->m_total = a_maxRects;
	a_parVars->m_minFill = a_minFill;
	for(int index=0; index < a_maxRects; ++index)
	{
		a_parVars->m_taken[index] = false;
		a_parVars->m_partition[index] = -1;
	}
}


RTREE_TEMPLATE
void RTREE_QUAL::PickSeeds(PartitionVars* a_parVars)
{
	int seed0=0, seed1=0;
	ELEMTYPEREAL worst, waste;
	ELEMTYPEREAL area[MAXNODES+1];

	for(int index=0; index<a_parVars->m_total; ++index)
	{
		area[index] = CalcRectVolume(&a_parVars->m_branchBuf[index].m_rect);
	}

	//worst = -a_parVars->m_coverSplitArea - 1;
	worst =-1e+300;//HY:large number need this init value
	for(int indexA=0; indexA < a_parVars->m_total-1; ++indexA)
	{
		for(int indexB = indexA+1; indexB < a_parVars->m_total; ++indexB)
		{
			Rect oneRect = CombineRect(&a_parVars->m_branchBuf[indexA].m_rect, &a_parVars->m_branchBuf[indexB].m_rect);
			waste = CalcRectVolume(&oneRect) - area[indexA] - area[indexB];
			if(waste > worst)
			{
				worst = waste;
				seed0 = indexA;
				seed1 = indexB;
			}
		}
	}
	Classify(seed0, 0, a_parVars);
	Classify(seed1, 1, a_parVars);
}


// Put a branch in one of the groups.
RTREE_TEMPLATE
void RTREE_QUAL::Classify(int a_index, int a_group, PartitionVars* a_parVars)
{
	ASSERT(a_parVars);
	ASSERT(!a_parVars->m_taken[a_index]);

	a_parVars->m_partition[a_index] = a_group;
	a_parVars->m_taken[a_index] = true;

	if (a_parVars->m_count[a_group] == 0)
	{
		a_parVars->m_cover[a_group] = a_parVars->m_branchBuf[a_index].m_rect;
	}
	else
	{
		a_parVars->m_cover[a_group] = CombineRect(&a_parVars->m_branchBuf[a_index].m_rect, &a_parVars->m_cover[a_group]);
	}
	a_parVars->m_area[a_group] = CalcRectVolume(&a_parVars->m_cover[a_group]);
	++a_parVars->m_count[a_group];
}


// Delete a data rectangle from an index structure.
// Pass in a pointer to a Rect, the tid of the record, ptr to ptr to root node.
// Returns 1 if record not found, 0 if success.
// RemoveRect provides for eliminating the root.
RTREE_TEMPLATE
bool RTREE_QUAL::RemoveRect(Rect* a_rect, const DATATYPE& a_id, Node** a_root)
{
	ASSERT(a_rect && a_root);
	ASSERT(*a_root);

	Node* tempNode;
	ListNode* reInsertList = NULL;

	if(!RemoveRectRec(a_rect, a_id, *a_root, &reInsertList))
	{
		// Found and deleted a data item
		// Reinsert any branches from eliminated nodes
		while(reInsertList)
		{
			tempNode = reInsertList->m_node;

			for(int index = 0; index < tempNode->m_count; ++index)
			{
				InsertRect(&(tempNode->m_branch[index].m_rect),
				           tempNode->m_branch[index].m_data,
				           a_root,
				           tempNode->m_level);
			}

			ListNode* remLNode = reInsertList;
			reInsertList = reInsertList->m_next;

			FreeNode(remLNode->m_node);
			FreeListNode(remLNode);
		}

		// Check for redundant root (not leaf, 1 child) and eliminate
		if((*a_root)->m_count == 1 && (*a_root)->IsInternalNode())
		{
			tempNode = (*a_root)->m_branch[0].m_child;

			ASSERT(tempNode);
			FreeNode(*a_root);
			*a_root = tempNode;
		}
		return false;
	}
	else
	{
		return true;
	}
}


// Delete a rectangle from non-root part of an index structure.
// Called by RemoveRect.  Descends tree recursively,
// merges branches on the way back up.
// Returns 1 if record not found, 0 if success.
RTREE_TEMPLATE
bool RTREE_QUAL::RemoveRectRec(Rect* a_rect, const DATATYPE& a_id, Node* a_node, ListNode** a_listNode)
{
	ASSERT(a_rect && a_node && a_listNode);
	ASSERT(a_node->m_level >= 0);

	if(a_node->IsInternalNode())   // not a leaf node
	{
		for(int index = 0; index < a_node->m_count; ++index)
		{
			if(Overlap(a_rect, &(a_node->m_branch[index].m_rect)))
			{
				if(!RemoveRectRec(a_rect, a_id, a_node->m_branch[index].m_child, a_listNode))
				{
					if(a_node->m_branch[index].m_child->m_count >= MINNODES)
					{
						// child removed, just resize parent rect
						a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
					}
					else
					{
						// child removed, not enough entries in node, eliminate node
						ReInsert(a_node->m_branch[index].m_child, a_listNode);
						DisconnectBranch(a_node, index); // Must return after this call as count has changed
					}
					return false;
				}
			}
		}
		return true;
	}
	else     // A leaf node
	{
		for(int index = 0; index < a_node->m_count; ++index)
		{
			if(a_node->m_branch[index].m_child == reinterpret_cast<Node*> (a_id))
			{
				DisconnectBranch(a_node, index); // Must return after this call as count has changed
				return false;
			}
		}
		return true;
	}
}


// Decide whether two rectangles overlap.
RTREE_TEMPLATE
bool RTREE_QUAL::Overlap(Rect* a_rectA, Rect* a_rectB)
{
	ASSERT(a_rectA && a_rectB);

	for(int index=0; index < NUMDIMS; ++index)
	{
		if (a_rectA->m_min[index] > a_rectB->m_max[index] ||
		        a_rectB->m_min[index] > a_rectA->m_max[index])
		{
			return false;
		}
	}
	return true;
}

// Decide whether a point is in a rectangle.
RTREE_TEMPLATE
bool RTREE_QUAL::Overlap(Point* a_point, Rect* a_rect)
{
	ASSERT(a_point && a_rect);

	for(int index=0; index < NUMDIMS; ++index)
	{
		if (a_point->m[index] > a_rect->m_max[index] ||
		        a_rect->m_min[index] > a_point->m[index])
		{
			return false;
		}
	}
	return true;
}


// Add a node to the reinsertion list.  All its branches will later
// be reinserted into the index structure.
RTREE_TEMPLATE
void RTREE_QUAL::ReInsert(Node* a_node, ListNode** a_listNode)
{
	ListNode* newListNode;

	newListNode = AllocListNode();
	newListNode->m_node = a_node;
	newListNode->m_next = *a_listNode;
	*a_listNode = newListNode;
}


// Search in an index tree or subtree for all data retangles that overlap the argument rectangle.
RTREE_TEMPLATE
bool RTREE_QUAL::Search(Node* a_node, Rect* a_rect, int& a_foundCount, void* a_context, bool/* __cdecl*/ a_resultCallback(DATATYPE a_data, void* a_context))
{
	ASSERT(a_node);
	ASSERT(a_node->m_level >= 0);
	ASSERT(a_rect);

	if(a_node->IsInternalNode())   // This is an internal node in the tree
	{
		for(int index=0; index < a_node->m_count; ++index)
		{
			if(Overlap(a_rect, &a_node->m_branch[index].m_rect))
			{
				if(!Search(a_node->m_branch[index].m_child, a_rect, a_foundCount, a_context, a_resultCallback))
				{
					return false; // Don't continue searching
				}
			}
		}
	}
	else     // This is a leaf node
	{
		for(int index=0; index < a_node->m_count; ++index)
		{
			if(Overlap(a_rect, &a_node->m_branch[index].m_rect))
			{
				DATATYPE& id = a_node->m_branch[index].m_data;

				// NOTE: There are different ways to return results.  Here's where to modify
				if(a_resultCallback!=NULL)
				{
					++a_foundCount;
					if(!a_resultCallback(id, a_context))
					{
						return false; // Don't continue searching
					}
				}
			}
		}
	}

	return true; // Continue searching
}



RTREE_TEMPLATE
bool RTREE_QUAL::Search(Node* a_node, Point* a_point, int& a_foundCount, void* a_context, bool/* __cdecl*/ a_resultCallback(DATATYPE a_data, void* a_context))
{
	ASSERT(a_node);
	ASSERT(a_node->m_level >= 0);
	ASSERT(a_point);

	if(a_node->IsInternalNode())   // This is an internal node in the tree
	{
		for(int index=0; index < a_node->m_count; ++index)
		{
			if(Overlap(a_point, &a_node->m_branch[index].m_rect))
			{
				if(!Search(a_node->m_branch[index].m_child, a_point, a_foundCount, a_context, a_resultCallback))
				{
					return false; // Don't continue searching
				}
			}
		}
	}
	else     // This is a leaf node
	{
		for(int index=0; index < a_node->m_count; ++index)
		{
			if(Overlap(a_point, &a_node->m_branch[index].m_rect))
			{
				DATATYPE& id = a_node->m_branch[index].m_data;

				// NOTE: There are different ways to return results.  Here's where to modify
				if(a_resultCallback!=NULL)
				{
					++a_foundCount;
					if(!a_resultCallback(id, a_context))
					{
						return false; // Don't continue searching
					}
				}
			}
		}
	}

	return true; // Continue searching
}

#define ABS(a)   (a)>0?(a):-(a)
RTREE_TEMPLATE
double RTREE_QUAL::maxd(const ELEMTYPE p[],const ELEMTYPE a_min[], const ELEMTYPE a_max[], int d_index)
{
	ELEMTYPE d_i_min=p[d_index]-a_min[d_index],d_i_max=p[d_index]-a_max[d_index];
	d_i_min=ABS(d_i_min);
	d_i_max=ABS(d_i_max);
	const ELEMTYPE *a_2d[2];
	a_2d[0]=a_min;
	a_2d[1]=a_max;
	int i_num;
	if(d_i_min>d_i_max)
		i_num=1;
	else
		i_num=0;
	ELEMTYPE sum=0,s_max,t;
	for(int i=0; i<NUMDIMS; i++)
	{
		s_max=0;
		if(i==d_index)
		{
			t=a_2d[i_num][i]-p[i];
			s_max=ABS(t);
		}
		else
			for(int j=0; j<2; j++)
			{
				t=a_2d[j][i]-p[i];
				t=ABS(t);
				if(t>s_max)
					s_max=t;
			}
		sum+=s_max*s_max;
	}
	return sqrt(sum);
}
RTREE_TEMPLATE
double RTREE_QUAL::cal_mindist(const ELEMTYPE p[],const ELEMTYPE a_min[],const ELEMTYPE a_max[])
{

	const ELEMTYPE *a_2d[2];
	a_2d[0]=a_min;
	a_2d[1]=a_max;
	ELEMTYPE sum=0,s_min,t;
	for(int i=0; i<NUMDIMS; i++)
	{
		s_min=std::numeric_limits<ELEMTYPE>::max();

		for(int j=0; j<2; j++)
		{
			t=a_2d[j][i]-p[i];
			t=ABS(t);
			if(t<s_min)
				s_min=t;
		}
		if (a_2d[0][i]<=p[i]&&a_2d[1][i]>=p[i])
			s_min=0;

		sum+=s_min*s_min;
	}
	return sqrt(sum);
}
RTREE_TEMPLATE
double RTREE_QUAL::cal_minmaxdist(const ELEMTYPE p[],const ELEMTYPE a_min[],const ELEMTYPE a_max[])
{
	ELEMTYPE min=std::numeric_limits<ELEMTYPE>::max();
	ELEMTYPE t;
	for(int i=0; i<NUMDIMS; i++)
	{

		t=maxd(p, a_min,a_max, i);
		if(t<min)
			min=t;
	}
	return min;
}
RTREE_TEMPLATE
double RTREE_QUAL::cal_mindist(const ELEMTYPE p[],const Rect* a_rect)
{

	const ELEMTYPE* a_min=a_rect->m_min;
	const ELEMTYPE* a_max=a_rect->m_max;
	return cal_mindist(p, a_min,a_max);
}
RTREE_TEMPLATE
double RTREE_QUAL::cal_minmaxdist(const ELEMTYPE p[],const Rect* a_rect)
{

	const ELEMTYPE* a_min=a_rect->m_min;
	const ELEMTYPE* a_max=a_rect->m_max;
	return cal_minmaxdist(p, a_min,a_max);
}


RTREE_TEMPLATE
int RTREE_QUAL::NNSearch(const ELEMTYPE a_p[])
{
	double minmaxdist,mindist;
	Branch* near_p=NULL;
	double near_d=std::numeric_limits<double>::max();
	std::multimap<double,tree_Node,std::less<double> > rbtree;
	typename std::multimap<double,tree_Node>::iterator pos;
	double mm=std::numeric_limits<double>::max();
	tree_Node temp_Node;
	temp_Node.m_mm=mm;
	temp_Node.m_node=m_root;
	rbtree.insert(std::make_pair(0,temp_Node));
	std::pair<double,tree_Node> temp;
	rbtree.size();
	while(!rbtree.empty())
	{
		temp=*(rbtree.begin());
		rbtree.erase(rbtree.begin());
		if(temp.second.m_node->IsInternalNode())   // This is an internal node in the tree
		{

			for(int index=0; index < temp.second.m_node->m_count; ++index)
			{
				minmaxdist=cal_minmaxdist(a_p,&temp.second.m_node->m_branch[index].m_rect);
				mindist=cal_mindist(a_p,&temp.second.m_node->m_branch[index].m_rect);
				if(minmaxdist<mm)
					mm=minmaxdist;
				temp_Node.m_mm=minmaxdist;
				temp_Node.m_node=temp.second.m_node->m_branch[index].m_child;
				rbtree.insert(std::make_pair(mindist,temp_Node));
			}
			for(pos = rbtree .begin(); pos!=rbtree.end() ; pos++ )
			{

				bool flag=rbtree.empty();
				while(!rbtree.empty()&&mm<pos->first)
				{

					pos=rbtree.erase(pos);

				}
				if(rbtree.empty())
					break;



			}
		}
		else
		{

			for(int index=0; index < temp.second.m_node->m_count; ++index)
			{
				double distance= cal_mindist(a_p,&temp.second.m_node->m_branch[index].m_rect);
				if(distance<near_d)
				{
					near_d=distance;
					near_p=&temp.second.m_node->m_branch[index];
					if(mm>near_d)
						mm=near_d;
				}

			}

		}
	}
	return near_p->m_data;
}

RTREE_TEMPLATE
inline typename RTREE_QUAL::Point RTREE_QUAL::Point::operator+(const typename RTREE_QUAL::Point& a)const
{
	typename RTREE_QUAL::Point t;
	for (int i=0; i<NUMDIMS; i++)
		t.m[i]=this->m[i]+a.m[i];
	return t;
};

RTREE_TEMPLATE
inline typename RTREE_QUAL::Point RTREE_QUAL::Point::operator-()const
{
	typename RTREE_QUAL::Point t;
	for (int i=0; i<NUMDIMS; i++)
		t.m[i]=-(this->m[i]);
	return t;
};
RTREE_TEMPLATE
inline typename RTREE_QUAL::Point RTREE_QUAL::Point::operator-(const typename RTREE_QUAL::Point& a)const
{
	return *this+(-a);
};
RTREE_TEMPLATE
inline typename RTREE_QUAL::Point RTREE_QUAL::Point::operator*(const double & a)const
{
	typename RTREE_QUAL::Point t;
	for (int i=0; i<NUMDIMS; i++)
		t.m[i]=a*(this->m[i]);
	return t;
};


#undef RTREE_TEMPLATE
#undef RTREE_QUAL

