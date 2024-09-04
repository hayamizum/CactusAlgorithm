#include <bits/stdc++.h>
#include <cfloat>
#define INT_MAX 2147483647
#define MAX_INT INT_MAX
using namespace std;

string file_name;
int mode; // mode 1 requires csv file, 2 just input.
double DistanceMatrix[1000][1000];
double Weight[1000][1000];
double SimilarityParameter=0.1;
bool Adjacent[1000][1000];
int n,n_Auxilliary=0;
vector<int>Now_Graph;

void CSV_READ() {
    ifstream file(file_name.c_str());
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
			
	//queue initialization
	for(int i=0;i<n;i++)
		Now_Graph.push_back(i);
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
    
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			Weight[i][j]=Weight[j][i]=DistanceMatrix[i][j];
}

void Compactification(int v,double cx){
	n_Auxilliary++;
	int new_vertex=n-1+n_Auxilliary;
	Now_Graph.push_back(new_vertex);
	
	DistanceMatrix[v][new_vertex]=DistanceMatrix[new_vertex][v]=cx;
	Weight[v][new_vertex]=Weight[new_vertex][v]=cx;
	
	for(int i=0;i<new_vertex;i++){
		
		if(i==v)continue;
		
		DistanceMatrix[i][new_vertex]=DistanceMatrix[new_vertex][i]=DistanceMatrix[i][v]-cx;
		Weight[new_vertex][i]=Weight[i][new_vertex]=INT_MAX;
		
		if(Weight[v][i]!=INT_MAX){
			
			Weight[new_vertex][i]=Weight[i][new_vertex]=Weight[v][i]-cx;
			Weight[v][i]=Weight[i][v]=INT_MAX;
			
		}
	}
}

void Check_Compactification(){
	int Size=Now_Graph.size();
	cout<<"Size->"<<Size<<'\n';
	for(int i=0;i<Size;i++){
		int vertex_i=Now_Graph[i];
		cout<<"vertex_i->"<<vertex_i<<'\n';
		double CurrentCx=INT_MAX;
		for(int j=0;j<Size;j++){
			int vertex_j=Now_Graph[j];
			if(i==j)continue;
			for(int k=j+1;k<Size;k++){
				if(i==k)continue;
				int vertex_k=Now_Graph[k];
				double Cx=DistanceMatrix[vertex_i][vertex_j]+DistanceMatrix[vertex_i][vertex_k]-DistanceMatrix[vertex_j][vertex_k];
				Cx/=2.0;
				if(Cx==0){
					goto early_ending;
				}
				CurrentCx=min(CurrentCx,Cx);
			}
		}
		cout<<"CurrentCx->"<<CurrentCx<<'\n';
		Compactification(vertex_i,CurrentCx);
		early_ending:
			while(0);
			//ending
	}
}

int get_degree(int v){
	int degree=0;
	for(int i=0;i<Now_Graph.size();i++){
		if(i==v)continue;
		if(Weight[i][v]!=MAX_INT)degree++;
	}
	return degree;
}

void Check_Edge(int x,int y){
	for(int i=0;i<Now_Graph.size();i++){
		int vertex=Now_Graph[i];
		if(vertex==x||vertex==y)continue;
		if(DistanceMatrix[x][vertex]+DistanceMatrix[vertex][y]==Weight[x][y]){
			if(DistanceMatrix[x][vertex]==0||DistanceMatrix[y][vertex]==0)continue;
			Weight[x][y]=Weight[y][x]=MAX_INT;
			puts("clear");
			cout<<"x->"<<x<<" y->"<<y<<" z->"<<vertex<<'\n';
		}
	}
}

void Clear_Graph_Edge(){
	for(int i=0;i<Now_Graph.size();i++){
		int vertex_i=Now_Graph[i];
		for(int j=i+1;j<Now_Graph.size();j++){
			int vertex_j=Now_Graph[j];
			if(Weight[vertex_i][vertex_j]!=MAX_INT){
				Check_Edge(vertex_i,vertex_j);
			}
		}
	}
}

void Clear_Graph_Vertex(){
	set<int>Graph;
	for(int i=0;i<Now_Graph.size();i++){
		int vertex=Now_Graph[i];
		if(get_degree(vertex)>2){
			Graph.insert(vertex);
		}
	}
	Now_Graph.clear();
	for(set<int>::iterator it=Graph.begin();it!=Graph.end();it++){
		Now_Graph.push_back(*it);
	}
}

void test_weight(){
	for(int i=0;i<n+n_Auxilliary;i++){
		for(int j=0;j<n+n_Auxilliary;j++){
			cout<<Weight[i][j]<<" ";
		}
		cout<<endl;
	}
}

void test_dis(){
	for(int i=0;i<n+n_Auxilliary;i++){
		for(int j=0;j<n+n_Auxilliary;j++){
			cout<<DistanceMatrix[i][j]<<" ";
		}
		cout<<endl;
	}
}

void Algorithm(){
	Check_Compactification();
	test_weight();
	Clear_Graph_Edge();
	Clear_Graph_Vertex();
}

void Print_Graph(){
	
	for(int i=0;i<n+n_Auxilliary;i++){
		for(int j=0;j<n+n_Auxilliary;j++){
			if(Weight[i][j]!=INT_MAX&&Weight[i][j]!=0)
				cout<<Weight[i][j]<<" ";
			else
				cout<<'x'<<" ";
		}
		cout<<endl;
	}
	
}

int main() {
	
    DistanceMatrixReadIn();
    AdjacentInitialize();
    // You can now use DistanceMatrix for further calculations
    FloydWarshall();
    Algorithm();
    Print_Graph();
    
    return 0;
}
