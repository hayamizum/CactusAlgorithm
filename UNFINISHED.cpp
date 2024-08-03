#include <bits/stdc++.h>
#include <cfloat>
using namespace std;

string file_name;
int mode; // mode 1 requires csv file, 2 just input.
double DistanceMatrix[1000][1000];
double Weight[1000][1000];
double SimilarityParameter=0.1;
bool Adjacent[1000][1000];
int n,n_Auxilliary;

void CSV_READ() {
    ifstream file(file_name);
    if (!file.is_open()) {
        cerr << "Failed to open file: " << file_name << endl;
        exit(1);
    }

    string line;
    // Read the first line to get the size of the matrix
    getline(file, line);
    istringstream ss(line);
    ss >> n;

    // Read the next n lines to get the distance matrix
    for (int i = 0; i < n; ++i) {
        getline(file, line);
        istringstream rowStream(line);
        for (int j = 0; j < n; ++j) {
            rowStream >> DistanceMatrix[i][j];
            // If CSV uses commas as separators, use rowStream with a custom delimiter
            // char delimiter;
            // rowStream >> DistanceMatrix[i][j] >> delimiter;
        }
    }

    file.close();
}

void DistanceMatrixReadIn(){
	cout << "MODE select, enter 1 for csv mode, 2 for manual mode\ninput here: ";
    cin >> mode;
    cout << endl;
    
	if (mode == 1) {
        cout << "Please enter the path of the csv file (absolute or relative)\ninput here: ";
        cin >> file_name;
        cout << endl;
        CSV_READ();
        cout << "CSV file has been read successfully! Start calculating...\n";
    } else if (mode == 2) {
        cout << "Please first enter the number (n) of points\ninput here: ";
        cin >> n;
        cout << endl;
        cout << "Please enter " << n << "x" << n << " numbers below\ninput here:\n";

        int COUNT = n * n;
        while (COUNT--) cin >> DistanceMatrix[COUNT / n][COUNT % n];

        cout << "Completed! Start calculating...\n";
    } else {
        cout << "This parameter is not supported, please try again";
        exit(0);
    }
} 

void AdjacentInitialize(){
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			Adjacent[i][j]=(i!=j);
}

void FloydWarshall(){
	for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (DistanceMatrix[i][k] < INT_MAX && DistanceMatrix[k][j] < INT_MAX && DistanceMatrix[i][k] + DistanceMatrix[k][j] < DistanceMatrix[i][j]) {
                    DistanceMatrix[i][j] = DistanceMatrix[i][k] + DistanceMatrix[k][j];
                }
            }
        }
    }
    
    for(int i=0;i<n;i++){
    	for(int j=0;j<n;j++){
    		Weight[i][j]=Weight[i][j];
		}
	}
    
}

void Compactification(int x,vector<int> vertex){
	
	double MIN=DBL_MAX;
	int y,z;
	
	for(int i_y=0;i_y<vertex.size();i_y++){
		_y=vertex[i_y];
		if(Adjacent[x][_y]==0)continue;
		for(int i_z=i_y+1;i_z<n;i_z++){
			_z=vertex[i_z];
			if(Adjacent[x][_z]==0)continue;
			if(MIN>Weight[x][_y]+Weight[x][_z]-DistanceMatrix[_y][_z])){
				MIN=Weight[x][_y]+Weight[x][_z]-DistanceMatrix[_y][_z];
				y=_y;
				z=_z;
			}
		}
	}
	MIN/=2.0;
	if(MIN<0)return;
	
	n_Auxilliary++;
	NewAuxilliaryNumber=n+n_Auxilliary-1;//warning
	vertex.push_back(NewAuxilliaryNumber);
	
	Adjacent[NewAuxilliaryNumber][x]=Adjacent[x][NewAuxilliaryNumber]=true;
	Adjacent[NewAuxilliaryNumber][y]=Adjacent[y][NewAuxilliaryNumber]=true;
	Adjacent[z][NewAuxilliaryNumber]=Adjacent[NewAuxilliaryNumber][z]=true;
	Adjacent[x][y]=Adjacent[y][x]=false;
	Adjacent[x][z]=Adjacent[z][x]=false;
	
	Weight[x][NewAuxilliaryNumber]=Weight[NewAuxilliaryNumber][x]=MIN;
	Weight[y][NewAuxilliaryNumber]=Weight[NewAuxilliaryNumber][y]=Weight[x][y]-MIN;
	Weight[z][NewAuxilliaryNumber]=Weight[NewAuxilliaryNumber][z]=Weight[x][z]-MIN;
	Weight[x][y]=Weight[y][x]=DBL_MAX;
	Weight[x][z]=Weight[z][x]=DBL_MAX;
	
	for(int i=0;i<NewAuxilliaryNumber;i++){
		if(i==x||i==y||i==z){
			DistanceMatrix[i][NewAuxilliaryNumber]=DistanceMatrix[NewAuxilliaryNumber][i]=Weight[NewAuxilliaryNumber][i];
			continue;
		}
		
		Max=DistanceMatrix[x][i]-DistanceMatrix[x][NewAuxilliaryNumber];
		Max=max(Max,DistanceMatrix[y][i]-DistanceMatrix[y][NewAuxilliaryNumber]);
		Max=max(Max,DistanceMatrix[z][i]-DistanceMatrix[z][NewAuxilliaryNumber]);
		DistanceMatrix[i][NewAuxilliaryNumber]=DistanceMatrix[NewAuxilliaryNumber][i]=Max;
	}
	DistanceMatrix[NewAuxilliaryNumber][NewAuxilliaryNumber]=0;
	
}

void Compaction(vector<int> vertex){
	for(int i=0;i<vertex.size();i++)
		Compactification(vertex[i]);
}

bool CheckRedundantEdge(x,y,vertex){
	for(int i=0;i<vertex.size();i++){
		
		if(vertex[i]==x||vertex[i]==y)continue;
		
		double SUM=DistanceMatrix[x][vertex[i]]+DistanceMatrix[vertex[i]][y];
		double proportion=SUM/Weight[x][y];
		
		if(proportion<(1+SimilarityParameter)) return true;
		else return false;
		
	}
}

vector<int> GetG_(vector<int> vertex){
	
	for(int i=0;i<vertex.size();i++){
		for(int j=i+1;j<vertex.size();j++){
			if(Adjacent[vertex[i]][vertex[j]]){
				if(CheckRedundantEdge(vertex[i],vertex[j],vertex)){
					Weight[vertex[i]][vertex[j]]=Weight[vertex[j]][vertex[i]]=DBL_MAX;
					Adjacent[vertex[i]][vertex[j]]=Adjacent[vertex[j]][vertex[i]]=0;
				}
			}
		}
	}
	
	vector<int> G_;
	
	for(int i=0;i<vertex.size();i++){
		
		int NumberOfNeighbors=0;
		for(int j=0;j<vertex.size();j++){
			if(i==j)continue;
			if(Adjacent[vertex[i]][vertex[j]])NumberOfNeighbors++;
		}
		
		if(NumberOfNeighbors>2){
			G_.push_back(vertex[i]);
		}
		
	}
	
	return G_;
}

void MainAlgorithmPart(){
	
	vector<int> VertexInG_previous(n),VertexInG_now;
	iota(VertexInG_previous.begin(),VertexInG_previous.end(),0);
	
	Step2:
		
		Compaction(VertexInG_previous);
		VertexInG_now=GetG_(VertexInG_previous);
		if(VertexInG_now.size()>1){
			VertexInG_previous=VertexInG_now;
			goto Step2;
		}
		
		
	Step3:
		VertexInG_now=VertexInG_previous;
		
		

}

int main() {
	
    DistanceMatrixReadIn();
    AdjacentInitialize();
    // You can now use DistanceMatrix for further calculations
    FloydWarshall();
    MainAlgorithmPart();
    
    return 0;
}
