#define rowLength		63
#define columnLength	63
#define totalWeeks		832



struct matrix
{
    short location[3969][3969];
};
//creates seaIce node which holds the deviation amount and its adjacent cells
struct seaIceNode
{
	float deviation;
	bool valid;
};

//creates a node for each cell for one week
struct weeklyCellData
{
	seaIceNode cells[rowLength][columnLength];
};

//combines all data
struct allData
{
	weeklyCellData fulldata[totalWeeks];
};

//node in a graph
struct graphEntry
{
	std::string	cellLoc;
	int rowNum;
	int columnNum;
	float xbar;
	float sxx;
	bool land;
	bool visited;
	graphEntry *edgeTo;
	
};

//graph structure
struct graph
{
	graphEntry location[rowLength][columnLength];
};



