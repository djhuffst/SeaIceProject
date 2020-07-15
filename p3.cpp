#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <thread>
#include "structs.h"
#define inf 30000

//global variables
allData dataSet;
int landCells;
graph thisGraph;
graph thisGraph1;
graph thisGraph2;
int degree[3969];
int connectionSizes[3969];
int conSize = 0;
float adjMatrix[3969][3969];
float adjMatrix1[3969][3969];
float adjMatrix2[3969][3969];

float kArray[3];
int kCount=0;



int min(int a, int b)
{
	if(a < b)
	{
		return a;
	}
	else
	{
		return b;
	}
}

void FloydWarshall(float adjMatrix[][3969],std::string threshold)
{
	std::cout<<"\r";
	std::cout<<"                                 ";
	std::flush(std::cout);

	int percent=0;
	int count = 0;

	for(int k = 0;k < 3969;k++)
	{
		for(int i = 0;i < 3969;i++)
		{
			for(int j = 0;j < 3969; j++)
			{
				count++;
				if(count % 312617511 == 0)
				{
					percent++;
					std::cout<<"\r";
					std::cout<<"Calculating shortest Paths: "<<percent<<"%";
					std::flush(std::cout);
				}

				adjMatrix[i][j] = min(adjMatrix[i][j],(adjMatrix[i][k] + adjMatrix[k][j]));

			}
		}
	}
	float count0 = 0;

	float sum1 = 0;

	
	for(int i = 0; i < 3969; i++)
	{
		for(int j = 0; j < 3969; j++)
		{
			if(adjMatrix[i][j] != inf)
			{
				sum1 = sum1 + adjMatrix[i][j];
				count0++;
			}

		}
	}
	std::cout<<"\r";
	std::cout<<"                                    ";
	std::flush(std::cout);
	std::cout<<"\r";
	std::cout<<"Characteristic path length for "<<threshold<<": "<<(sum1 / count0)<<std::endl;
}

//sets matrix initial values to infinity
void clearMatrix()
{
	for(int i = 0; i < 3969; i++)
	{
		for(int j = 0; j < 3969; j++)
		{
			adjMatrix[i][j] = inf;
			adjMatrix1[i][j] = inf;
			adjMatrix2[i][j] = inf;
		}
	}
}

//checks adjacencies for connection
bool connected(int row,int column, int row2, int column2,graph * graphy)
{
	graphEntry *node = &graphy->location[row2][column2];
	while(node->edgeTo != NULL)
	{
		node = node->edgeTo;
		if(node->rowNum == row && node->columnNum == column)
		{
			return true;
		}
	}
	return false;
}

//calculates the clustering coefficient for a single vertex
float clustCoef(int row, int column,graph *graphy)
{
	float coef = 0;
	float vertices = 0;
	float edges = 0;
	graphEntry * node = new graphEntry();
	graphEntry * node2 = new graphEntry();
	
	node = &graphy->location[row][column];
	while(node->edgeTo != NULL)
	{ 
		vertices++;
		node = node->edgeTo;
	}

	node = &graphy->location[row][column];

	if(node->edgeTo != NULL)
	{
		node = node->edgeTo;
	}

	while(node->edgeTo != NULL)
	{
		node2 = node->edgeTo;
		if(connected(node->rowNum,node->columnNum,node2->rowNum,node2->columnNum,graphy))
		{
			edges++;
		}
		while(node2->edgeTo != NULL)
		{
			node2 = node2->edgeTo;
			if(connected(node->rowNum,node->columnNum,node2->rowNum,node2->columnNum,graphy))
			{
				edges++;
			}
			
		}
		node = node->edgeTo;
	}

	if(vertices > 1)
	{
		coef = ( 2 * edges) / (vertices * (vertices - 1));
	}
	else
	{
		coef = 0.0;
	}
	return coef;
}

//calculates the average Clustering Coefficient for a graph
float aveCoef(graph* graphy)
{
	float sum = 0.0;
	int row = 0;
	int column = 0;
	float count = 0.0;

	for(int j = 0;j < 3969;j++)
	{
		if(column == 63)
		{
			column = 0;
			row++;
		}
		if(dataSet.fulldata[0].cells[row][column].valid)
		{
			float cluster = clustCoef(row,column,graphy);
			count += 1;
			sum = sum + cluster;	
		}
		column++;
	}
	return (sum / count);
}

//Depth-Fist-Search of graph
void DFS(int row, int column, graph* graphy)
{
	conSize++;
	graphy->location[row][column].visited = true;

	graphEntry * tmpNode = new graphEntry();
	tmpNode = &graphy->location[row][column];

	while(tmpNode->edgeTo != NULL)
	{
		tmpNode = tmpNode->edgeTo;
		if(graphy->location[tmpNode->rowNum][tmpNode->columnNum].visited == false)
		{	
			DFS(tmpNode->rowNum,tmpNode->columnNum,graphy);
		}
	}
}


//traverses graph with DFS and counts connections
int conComponents(graph *graphy)
{
	int count = 0;
	int otherCount = 0;
	int column = 0;
	int row = 0;

	for(int k = 0;k < 3969;k++)
	{
		connectionSizes[k] = 0;
	}

	for(int i = 0;i < 3969;i++)
	{
		conSize = 0;
		
		if(column == 63)
		{
			column = 0;
			row++;
		}

		if(graphy->location[row][column].visited == false && !(graphy->location[row][column].land))
		{
			DFS(row,column,graphy);
			count++;
		}

		if(conSize > 0)
		{
			connectionSizes[conSize]++;
		}
		column++;
	}
	return count;	
}

//prints out the degree distribution
void degreeDist(graph* graphy,std::string threshold)
{

	std::cout<<""<<std::endl;
	std::cout<<""<<std::endl;
	std::cout<<"*******************************"<<std::endl;
	std::cout<<"	"<<threshold<<" threshold		   "<<std::endl;
	std::cout<<"*******************************"<<std::endl;
	int sNodeCounter = 0;
	int row = 0;
	int column = 0;
	std::string cell;

	float degreeCount = 0;
	float sum = 0;

	for(int j = 0;j < 3969;j++)
	{
		degree[j] = 0;
	}
	for(int i = 0;i < 3969;i++)
	{
		int degreeCount = 0;
		if(column == 63)
		{
			column = 0;
			row++;
		}

		graphEntry * thisNode = &graphy->location[row][column];
		while(thisNode->edgeTo != NULL)
		{
			degreeCount++;
			
			thisNode = thisNode->edgeTo;
		}
			if(thisNode->land == false)
			{
				int size = degree[degreeCount];
				degree[degreeCount] = size+1;
			}

		column++;
	}
	for(int j = 0;j < 3969;j++)
	{
		if(degree[j] != 0)
		{
			
			std::cout<<"Degree "<<j<<": ";
			for(int k = 0;k < degree[j];k += 10)
			{
				std::cout<<"*";
			}
			std::cout<<"("<<degree[j]<<")"<<std::endl;
			
			//random graph calculations
			degreeCount+=degree[j];
			sum += (degree[j] * j);
		}
	}
	kArray[kCount] = sum / degreeCount;
	kCount++;
}


//Pearson correlation coefficient
float correlation(int xPos,int yPos,int xPos2, int yPos2)
{
	
	float sxy = 0;
	float r = 0;
	
	//pull values from graph nodes
	float xbar=thisGraph.location[xPos][yPos].xbar;
	float ybar=thisGraph.location[xPos2][yPos2].xbar;
	float sxx=thisGraph.location[xPos][yPos].sxx;
	float syy=thisGraph.location[xPos2][yPos2].sxx;

	//checks too see if either value is a land node
	for(int k = 0;k < 832;k++)
	{
		sxy = sxy + ((dataSet.fulldata[k].cells[xPos][yPos].deviation - xbar) * (dataSet.fulldata[k].cells[xPos2][yPos2].deviation - ybar));
	}
	//calculate r value
	r = fabs(sxy / (sqrt(sxx * syy)));
	
	return r;

}

void createAdjacencies()
{
	std::cout<<"\r";
	std::cout<<"                                 ";
	std::flush(std::cout);

	int row1 = 0;
	int column1 = 0;
	int row2 = 0;
	int column2 = 0;
	int percentCount;
	int percentComplete = 0;

	for(int i = 0;i < 3969;i++)
	{	
		
		if(column1 == 63)
		{
			column1 = 0;
			row1++;
		}

		for(int j = 0;j < 3969;j++)
		{
			
			percentCount++;

			if(percentCount % 157529 == 0)
			{
				std::cout<<"\r";
				percentComplete++;
				std::cout<<"Building Graphs: "<<percentComplete<<"%";
				std::flush(std::cout);
			}
		
			if(column2 == 63)
			{
				column2 = 0;
				row2++;
			}
			if(j == 0)
			{
				row2 = 0;
			}
			if(!(row1 == row2 && column1 == column2) && dataSet.fulldata[0].cells[row1][column1].valid && dataSet.fulldata[0].cells[row2][column2].valid)
			{
				double r = correlation(row1,column1,row2,column2);

				//90% correlation
				if(r >= .90)
				{
					//creates an edge in adjacency matrix
					adjMatrix[(row1 * 63) + (column1) ][(row2 * 63) + (column2)] = 1;
					
					//adds adjacency node to adjacency list
					graphEntry *newNode = new graphEntry();
					newNode->rowNum = row2;
					newNode->columnNum = column2;
					newNode->edgeTo = NULL;

					graphEntry *edge=new graphEntry();
					edge=&thisGraph.location[row1][column1];

					while(edge->edgeTo!=NULL)
					{
						edge=edge->edgeTo;	
					}
					edge->edgeTo=newNode;
					
				}

				//92.5% correlation
				if(r>=.925)
				{
					adjMatrix1[(row1 * 63) + (column1) ][(row2 * 63) + (column2)] = 1;

					graphEntry *newNode=new graphEntry();
					newNode->rowNum=row2;
					newNode->columnNum=column2;
					newNode->edgeTo=NULL;

					graphEntry *edge=new graphEntry();
					edge=&thisGraph1.location[row1][column1];

					while(edge->edgeTo!=NULL)
					{
						edge=edge->edgeTo;	
					}
					edge->edgeTo=newNode;
					
				}
				
				//95% correlation
				if(r>=.95)
				{
					adjMatrix2[(row1 * 63) + (column1) ][(row2 * 63) + (column2)] = 1;

					graphEntry *newNode=new graphEntry();
					newNode->rowNum=row2;
					newNode->columnNum=column2;
					newNode->edgeTo=NULL;

					graphEntry *edge=new graphEntry();
					edge=&thisGraph2.location[row1][column1];

					while(edge->edgeTo!=NULL)
					{
						edge=edge->edgeTo;	
					}
					edge->edgeTo=newNode;
				}
			}
			column2++;
		}
	
		column1++;

	}
}


//creates and sets default values for graph
void createGraph(graph * graphy)
{
	std::cout<<"\r";
	std::cout<<"Building empty graphs";
	std::flush(std::cout);

	int row1=0;
	int column1=0;

	for(int i=0;i<3969;i++)
	{	
		if(column1==63)
		{
			column1=0;
			row1++;
		}
		graphEntry *node=new graphEntry();
		node->rowNum=row1;
		node->columnNum=column1;
		node->xbar=0;
		node->sxx=0;
		node->visited=false;
		node->edgeTo=NULL;
		graphy->location[row1][column1]=*node;
		column1++;
	}	
}


//calculate xbar and sxx for all nodes in the graph
void calcMeans()
{
	std::cout<<"\r";
	std::cout<<"Calculating xbar and sxx: ";
	std::flush(std::cout);

	int row=0;
	int column=0;
	bool isLand=true;

	for(int i=0;i<3969;i++)
	{
		float sum=0;
		float mean=0;
		float sxx=0;
		bool isLand=true;
		if(column==63)
		{
			column=0;
			row++;
		}
		if(dataSet.fulldata[0].cells[row][column].valid)
		{
			isLand=false;
			for(int j=0;j<832;j++)
			{
				sum=sum+dataSet.fulldata[j].cells[row][column].deviation;
			}
			mean=sum/832;

			for(int l=0;l<832;l++)
			{
				sxx=sxx+std::pow((dataSet.fulldata[l].cells[row][column].deviation-mean),2);	
			}
		}
		else
		{
			landCells++;
		}
		
		thisGraph.location[row][column].land=isLand;
		thisGraph.location[row][column].xbar=mean;
		thisGraph.location[row][column].sxx=sxx;

		thisGraph1.location[row][column].land=isLand;
		thisGraph1.location[row][column].xbar=mean;
		thisGraph1.location[row][column].sxx=sxx;

		thisGraph2.location[row][column].land=isLand;
		thisGraph2.location[row][column].xbar=mean;
		thisGraph2.location[row][column].sxx=sxx;	
		column++;
	}
	
}



//stores entire sea ice data
void storeData()
{
	std::cout<<"\r";
	std::cout<<"Storing data";
	std::flush(std::cout);

	weeklyCellData *weekly=new weeklyCellData();
	int weekCount=0;
	for(int i=1990;i<2006;i++)
	{	
		std::string year=std::to_string(i);
		for(int j=1;j<53;j++)
		{
			int rowCount=0;
			int columnCount=0;
			std::string week;
			if(j<10)
			{
				week="0"+std::to_string(j);
			}
			else
			{
				week=std::to_string(j);
			}
			std::ifstream inputFile( year+"/Beaufort_Sea_diffw"+week+"y"+year+"+landmask", std::ios::in | std::ios::binary );
			
			
			float dataIn = 0;

			for(int k=0;k<3969;k++)
			{
				if(columnCount==63)
				{
			
					columnCount=0;
					rowCount++;
				}
				
				inputFile.read( (char*)&dataIn, 4 ); 
		
				seaIceNode *cell=new seaIceNode();
				if(dataIn==168 || dataIn==157)
				{
					
					cell->valid=false;
					cell->deviation=dataIn;
					weekly->cells[rowCount][columnCount]=*cell;
				}
				else
				{
					
					cell->valid=true;
					cell->deviation=dataIn;
					weekly->cells[rowCount][columnCount]=*cell;
				}
				columnCount++;
			}
			
			dataSet.fulldata[weekCount]=*weekly;

			weekCount++;

			inputFile.close();
		}
	}
}



int main()
{

	storeData();
	createGraph(&thisGraph);
	createGraph(&thisGraph1);
	createGraph(&thisGraph2);
	clearMatrix();
	calcMeans();
	
	createAdjacencies();

	degreeDist(&thisGraph,".90");
	std::cout<<""<<std::endl;
	std::cout<<""<<std::endl;
	std::cout<<"Number of connected components at .90 threshold: "<<conComponents(&thisGraph)<<std::endl;
	std::cout<<"------------------------------------------------"<<std::endl;
	for(int j=0;j<3969;j++)
	{
		if(connectionSizes[j]!=0)
		{
			std::cout<<"Components with size "<<j<<": "<<connectionSizes[j]<<std::endl;
		}
	}

	degreeDist(&thisGraph1,".925");
	std::cout<<""<<std::endl;
	std::cout<<""<<std::endl;
	std::cout<<"Number of connected components at .925 threshold: "<<conComponents(&thisGraph1)<<std::endl;
	std::cout<<"-------------------------------------------------"<<std::endl;
	for(int j=0;j<3969;j++)
	{
		if(connectionSizes[j]!=0)
		{
			std::cout<<"Components with size "<<j<<": "<<connectionSizes[j]<<std::endl;
		}
	}
	
	degreeDist(&thisGraph2,".95");
	std::cout<<""<<std::endl;
	std::cout<<""<<std::endl;
	std::cout<<"Number of connected components at .95 threshold: "<<conComponents(&thisGraph2)<<std::endl;
	std::cout<<"------------------------------------------------"<<std::endl;
	for(int j=0;j<3969;j++)
	{
		if(connectionSizes[j]!=0)
		{
			std::cout<<"Components with size "<<j<<": "<<connectionSizes[j]<<std::endl;
		}
	}

	std::cout<<""<<std::endl;
	std::cout<<"-------------- Derived Graphs --------------"<<std::endl;
	std::cout<<""<<std::endl;
	std::cout<<"Average clustering coefficient .9: "<<aveCoef(&thisGraph)<<std::endl;
	std::cout<<"Average clustering coefficient .925: "<<aveCoef(&thisGraph1)<<std::endl;
	std::cout<<"Average clustering coefficient .95: "<<aveCoef(&thisGraph2)<<std::endl;
	std::cout<<""<<std::endl;

	std::thread thread3 {FloydWarshall, adjMatrix, ".9"};
	std::thread thread4 {FloydWarshall, adjMatrix1, ".925"};
	std::thread thread5 {FloydWarshall, adjMatrix2, ".95"};

	thread5.join();
	thread4.join();
	thread3.join();

	std::cout<<""<<std::endl;
	std::cout<<"--------------- Random Graph ---------------"<<std::endl;
	std::cout<<"Random Clustering coefficient .9: "<<kArray[0]/3186<<std::endl;
	std::cout<<"Random Clustering coefficient .925: "<<kArray[1]/3186<<std::endl;
	std::cout<<"Random Clustering coefficient .95: "<<kArray[2]/3186<<std::endl;

	std::cout<<""<<std::endl;
	std::cout<<"Random Characteristic Path Length .9: "<<log2(3186)/log2(kArray[0])<<std::endl;
	std::cout<<"Random Characteristic Path Length .925: "<<log2(3186)/log2(kArray[1])<<std::endl;
	std::cout<<"Random Characteristic Path Length .95: "<<log2(3186)/log2(kArray[2])<<std::endl;

 	return 0;
}

