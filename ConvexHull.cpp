#include <iostream>
#include <bits/stdc++.h>
#include <conio.h>

using namespace std;

/// A basic structure defining the coordinates of a point in the cartesian plane.
/*!
<pre>
 It also stores the index of the points, the sequence in which they were entered by the user.
 It defines the less than < operator to use for sort() function.
</pre>
*/
struct Point{
  /*! Stores the x-coordinate of a point. */
  int x;
  /*! Stores the y-coordinate of a point. */
  int y;
  /*! Stores the index of the point as entered by user. */
  int index;
  /*! Defining < Operator for Point structure. */
  bool operator < (const Point &p) const{
		if(x<p.x){
            return true;
		}
		else if(x==p.x&&y<=p.y){
            return true;
		}
		return false;
	}
};

/// A global point used for sorting points w.r.t. the first point.
/// Used for compare() function in qsort() function.
Point pref;

/// Class to generate Convex Hull of a set of points.
/*!
<pre>
It contains different functions to generate Convex Hull.
of a set of points using 3 different algorithms:
Jarvis's March
Graham's Scan
Andrew's Algorithm
It also contains utility functions to be used to implement these algorithms.
</pre>*/
class ConvexHull{
public:

/// Constructor
/*!
<pre>
 Constructor added in case needed in future.
</pre>
*/
ConvexHull(){}

/// Destructor
/*!
<pre>
 Destructor added in case needed in future.
</pre>
*/
~ConvexHull(){}



/// A function to create .ch file to be used by Visualizer.
/*!
<pre>
 The function takes a set of points and their convex hull as input and generates a .ch file
 in the required format to be used by Visualizer.
</pre>
\param vertices an array of point structures.
\param hull a vector<Point> containing elements of the Convex Hull of vertices.
\param n an integer. Size of vertices array.
\sa Point
*/
void toVisualizer(Point *vertices, vector<Point> hull,int n){
    int i;
    fstream fin;
    fin.open("hull.ch",ios::out);
    fin<<"CH"<<endl;
    fin<<n<<" "<<hull.size()<<endl;
    for(i=0;i<n;i++){
        fin<<vertices[i].x<<" "<<vertices[i].y<<" "<<"0"<<endl;
    }
    for(i=0;i<hull.size();i++){
        fin<<hull[i].index<<" ";
    }
    fin.close();
}



/// A function to find the orientation of the ordered triplet (p0, p1, p3).
/*!
<pre>
 It checks whether p0 -> p1 -> p2 makes a counter-clockwise or clockwise turn.
 It uses the crossproduct (p1-p0)X(p2-p0) to check the direction in which it turns.
</pre>
\param p0 a Point structure.
\param p1 a Point structure.
\param p2 a Point structure.
\return -1 if p0 -> p1 -> p2 makes a counter-clockwise turn.
\return 1 if p0 -> p1 -> p2 makes a clockwise turn.
\return 0 if the three points are collinear.
\sa Point, compare(), JarvisMarch(), GrahamScan(), AndrewAlgo()
*/
static int direction(Point p0, Point p1, Point p2){
    int crossproduct;

    crossproduct=(((p1.x-p0.x)*(p2.y-p0.y))-((p2.x-p0.x)*(p1.y-p0.y)));

    if(crossproduct==0){return 0;}

    else if(crossproduct<0){return 1;}

    else {return -1;}
}



/// A function to swap two points (entries), mainly in a Point array.
/*!
<pre>
 Used mainly for implementing sorting algorithms for a Point array.
</pre>
\param p1 a Point structure reference.
\param p2 a Point structure reference.
\sa Point, GrahamScan()
*/
void swapPoint(Point &p1, Point &p2){
    Point temp;
    temp=p1;
    p1=p2;
    p2=temp;
}



/// A function to find and return the square of the distance between two Points.
/*!
<pre>
 Used to decide which point to keep in case two points have same polar angles.
</pre>
\param p1 a Point structure.
\param p2 a Point structure.
\return dist - an integer value. Square of distance between p1 and p2.
\sa Point, compare(), JarvisMarch()
*/
static int dist(Point p1, Point p2){
    int dist;
    dist=(p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y);
    return dist;
}



/// A function to return the point which is next to the topmost point (second from top) in a stack.
/*!
<pre>
 Used in Graham's Scan algorithm.
</pre>
\param stck a reference to a vector<Point>, a vector of Point structures (type).
\return previous - a Point structure. The element next to the top element in stck.
\sa Point, GrahamScan()
*/
Point prev(vector<Point> &stck){
    Point current=stck.back();
    stck.pop_back();
    Point previous=stck.back();
    stck.push_back(current);
    return previous;
}



/// A function used by library function qsort to compare two Point structures' polar order with respect to the bottom-most point.
/*!
<pre>
 It also uses distance function, in case two points have same polar order.
 In such cases the nearest point is kept first.
</pre>
\param p1 a pointer to a Point structure.
\param p2 a pointer to a Point structure.
\return value - an integer
\sa Point, GrahamScan(), direction(), dist()
*/
static int compare(const void *p1, const void *p2){
    int value;
    Point *temp1=(Point *)p1;
    Point *temp2=(Point *)p2;

    // Find orientation
    int o=direction(pref, *temp1, *temp2);
    if(o==0){
        value=(dist(pref,*temp2)>=dist(pref,*temp1))?-1:1;
    }
    else{
        value=(o==-1)?-1: 1;
    }
    return value;
}




/// Function to return convex hull of a set of n points using Jarvis's March algorithm.
/*!
<pre>
<b>PSEUDOCODE:</b>
JarvisMarch(S)                                             // S is set of points whose Convex Hull is to be generated.
   leftmost = leftmost point in S
   i = 0
   current = leftmost
   repeat:
      P[i] = current                                       // P is set of points in Convex Hull.
      next = current+1 in set                              // initial endpoint for a candidate edge on the hull.
      for j from 0 to |S|:
         if (S[j] is on right of line from P[i] to next)   // find out point such that all triplets (current,S[j],next) are clockwise.
            next = S[j]                                    // found greater left turn, update endpoint.
      i = i+1
      current = next
   until current == leftmost                                // wrapped around to first hull point.
   return P

<b>TIME COMPLEXITY:</b>
Let n be number of total points entered and
    h be number of points in the Convex Hull.
The inner loop checks every point in the set S.
The outer loop repeats for each point on the hull.
Hence the total run time is O(n*h).
Worst case time complexity is O(n^2).
</pre>
\param vertices an array of Point structures.
\param n an integer value. Size of vertices array.
\return hullPoints - a vector<Point> (vector of Point structures). Contains the set of points in the Convex Hull of vertices.
\sa Point, direction(), dist()
*/
vector<Point> JarvisMarch(Point vertices[],int n){
    // minimum three points needed.
    if(n<3)return vector<Point> ();

    int i,current,next,leftmost=0,hullcount=0;

    // initialize result hull vector.
    vector<Point> hullPoints;

    // find leftmost point.
    for(i =0;i<n;i++){
        if(vertices[i].x<vertices[leftmost].x){
            leftmost=i;
        }
        if(vertices[i].x==vertices[leftmost].x){
            leftmost=(vertices[i].y<vertices[leftmost].y)?i:leftmost;
        }
    }

    current=leftmost;

    // starting from leftmost point, keep moving finding the most counter-clockwise point.
    // this loop runs O(h) times where h is number of points in result or output.
    // runs till we don't come to the first point.
    while(current!=leftmost||hullcount==0){
        // add current point to the result hull, and increase count of hull points.
        hullPoints.push_back(vertices[current]);
        hullcount++;

        // search for a point next, such that the direction of all triplets (current,vertices[i],next)
        // is counter-clockwise for all points i.
        // if any point i is more counter-clockwise, then update next to i;
        next=(current+1)%n;
        for(i=0;i<n;i++){
            // if i is more counter-clockwise then update next.
            if(direction(vertices[current],vertices[i],vertices[next])==-1){
                next=i;
            }
            if(direction(vertices[current],vertices[i],vertices[next])==0){
                next=dist(vertices[current],vertices[i])<=dist(vertices[current],vertices[next])?next:i;
            }
        }
        // next is the most counter-clockwise point now w.r.t. current.
        // set current to next for the next iteration, so next is added to the result hull.
        current=next;
    }

    /// Returns the resultant Convex Hull.
    /// Points in result will be listed in counter-clockwise order.
    return hullPoints;
}



/// Function to return convex hull of a set of n points using Graham's Scan algorithm.
/*!
<pre>
<b>PSEUDOCODE:</b>
GrahamScan(S)                                                // S is set of points whose Convex Hull is to be generated.
  first = lowermost point in S
  swap S[0] with first
  sort S by polar angle with S[0]
  if (two points have same polar angles)
    keep only the one farthest from S[0], save in M
  P[0] = S[0]                                                // P is the stack which has points in the Convex Hull.
  P[1] = S[1]                                                // put the first three points of S in stack P.
  P[2] = S[2]
  for i=3 to |M|:                                            // process the remaining |M| points to find next valid elements.
      while points (next-to-top in P,top of P,M[i])
            are not counter-clockwise:
          pop one element from stack P                       // keep removing elements from stack till next-to-top,top and M[i] make a non-left turn.
      push M[i] to P                                         // put the next valid element to stack.
  return P

<b>TIME COMPLEXITY:</b>
Let n be number of total points entered.
The algorithm takes O(nLogn) time if we use a O(nLogn) sorting algorithm.
The first step (finding the bottom-most point) takes O(n) time.
The second step (sorting points) takes O(nLogn) time.
The third step (keeping only points with different polar angles) takes O(n) time.
The step to process points one by one takes O(n) time, assuming that the stack operations take O(1) time.
Overall complexity is O(n) + O(nLogn) + O(n) + O(n) which is O(nLogn)
</pre>
\param vertices an array of Point structures.
\param n an integer value. Size of vertices array.
\return hull - a vector<Point> (vector of Point structures). Contains the set of points in the Convex Hull of vertices.
\sa Point, swapPoint(), compare(), prev(), direction()
*/
vector<Point> GrahamScan(Point vert[],int n){
    // minimum three points needed.
    if(n<3)return vector<Point> ();

    int i,first=0,modsize=1;
    Point vertices[n],mod[n];
    for(i=0;i<n;i++){
        vertices[i]=vert[i];
    }

    // initialize result hull stack.
    vector<Point> hull;

    // find the bottom-most point.
    // in case of a tie, take the leftmost point.
    for(i=0;i<n;i++){
        if(vertices[i].y<vertices[first].y){
            first=i;
        }
        if(vertices[i].y==vertices[first].y){
            first=(vertices[i].x<vertices[first].x)?i:first;
        }
    }

    // swap the bottom-most point with the first point in set
    swapPoint(vertices[0],vertices[first]);

    // sort points according to polar order w.r.t. first point
    pref=vertices[0];
    qsort(&vertices[1],n-1,sizeof(Point),compare);

    mod[0]=vertices[0];

    // make new set of points with unique polar angles.
    // in case of tie, keep point with the largest distance from first point.
    for(i=1;i<n;i++){
        while((i<n-1)&&(direction(vertices[0],vertices[i],vertices[i+1])==0)){
            i++;
        }
        mod[modsize]=vertices[i];
        modsize++;
    }

    // if modified set of points has less than 3 points, convex hull is not possible.
    if(modsize<3)return vector<Point> ();

    // push first three points in empty stack.
    hull.push_back(mod[0]);
    hull.push_back(mod[1]);
    hull.push_back(mod[2]);

    // processing remaining modsize-3 points.
    for(i=3;i<modsize;i++){
        while(direction(prev(hull),hull.back(),mod[i])!=-1){
            // keep removing top point (of stack) if
            //next-to-top,top,mod[i] are counter-clockwise oriented.
            hull.pop_back();
        }
        // add valid point (next point in convex hull) to the stack.
        hull.push_back(mod[i]);
    }

    /// Returns the resultant Convex Hull.
    /// Points in result will be listed in clockwise order (from stack-top to bottom).
    return hull;
}



/// Function to return convex hull of a set of n points using Andrew's algorithm.
/*!
<pre>
<b>PSEUDOCODE:</b>
AndrewAlgo(S)                                                // S is set of points whose Convex Hull is to be generated.
  sort S by x coordinate
  for i = 1 to n:                                            // P is set of points in Convex Hull.
    while P contains at least two points and the sequence of last two points
            of P and the point S[i] does not make a counter-clockwise turn:
        remove the last point from P
    append S[i] to P                                         // put the next valid element in convex hull.
  lowersize = size of elements in P
  for i = n to 1:
    while P contains at least lowersize+1 points and the sequence of last two points
            of P and the point S[i] does not make a counter-clockwise turn:
        remove the last point from P
    append S[i] to P
  return P

<b>TIME COMPLEXITY:</b>
Let n be number of total points entered.
The algorithm takes O(nLogn) time if we use a O(nLogn) sorting algorithm.
The first step (sorting points) takes O(nLogn) time.
The second and third step (constructing lower and upper hulls of the points) both take O(n) time.
Overall complexity is O(nLogn) + O(n) + O(n) which is O(nLogn).
</pre>
\param vert an array of Point structures.
\param n an integer value. Size of vert array.
\return hull - a vector<Point> (vector of Point structures). Contains the set of points in the Convex Hull of vertices.
\sa Point, direction()
*/
vector<Point> AndrewAlgo(Point vert[],int n){
    int i,hsize=0;
    vector<Point> vertices;
    for(i=0;i<n;i++){
        vertices.push_back(vert[i]);
    }

    // initialize result hull vector.
    vector<Point> hull;

    // sort points according to x coordinate.
    sort(vertices.begin(),vertices.end());


    // build the lower part of the Convex Hull.
    for(i=0;i<n;i++){
    // check if the next point in set makes a clockwise turn.
    // if yes then keep removing the last point from result hull.
        while((hsize>=2)&&(direction(hull.at(hsize-2),hull.at(hsize-1),vertices.at(i))!=-1)){
            hull.pop_back();
            hsize--;
		}
        // add the next valid point in result hull.
        hull.push_back(vertices[i]);
        hsize++;
	}

	int lowersize=hsize;

    // build the upper part of the Convex Hull.
	for(i=n-2;i>=0;i--){
    // check if the next point in set makes a clockwise turn.
    // if yes then keep removing the last point from result hull.
        while((hsize>=lowersize+1)&&(direction(hull.at(hsize-2),hull.at(hsize-1),vertices.at(i))!=-1)){
            hull.pop_back();
            hsize--;
		}
        // add the next valid point in result hull.
		hull.push_back(vertices[i]);
		hsize++;
	}

    // resize as the first point is also the last, hence redundant.
    hull.resize(hsize-1);

    /// Returns the resultant Convex Hull.
    /// Points in result will be listed in counter-clockwise order.
    return hull;
}
};

int main()
{
    int n,i,choice=0;
    cin>>n;
    ConvexHull obj;
    vector<Point> Hull;
    x:
    Point points[n];
    cout<<"Enter the points with index starting from 0."<<endl;
    for(i=0;i<n;i++){
        cin>>points[i].index>>points[i].x>>points[i].y;
    }

    //Point points[] = {{0, 3, 0}, {2, 2, 1}, {3, 3, 2}, {3, 2, 3},
    //               {3, 0, 4}, {1, 1, 5}, {3, 1, 6}, {2, 1, 7}, {0, 0, 8}};

    /*cout<<"Enter 1-4 to choose which algorithm to use."<<endl;
    cout<<"1. Jarvis's March."<<endl;
    cout<<"2. Graham's Scan."<<endl;
    cout<<"3. Andrew's Algorithm."<<endl;
    cout<<"4. Exit."<<endl;
    cin>>choice;*/

    //while(choice!=6&&choice!=4){

       // if(choice==1){
            Hull=obj.JarvisMarch(points, n);
            cout<<"JARVIS'S MARCH RESULTANT CONVEX HULL :\n";
            for(i=0;i<Hull.size();i++)
                cout<<Hull[i].index<<" ("<<Hull[i].x<<", "<<Hull[i].y<<")"<<endl;
        //}


        /*if(choice==2){
            Hull=obj.GrahamScan(points,n);
            cout<<"GRAHAM'S SCAN RESULTANT CONVEX HULL :\n";
            for(i=0;i<Hull.size();i++)
                cout<<Hull[i].index<<" ("<<Hull[i].x<<", "<<Hull[i].y<<")"<<endl;
        }


        if(choice==3){
            Hull=obj.AndrewAlgo(points,n);
            cout<<"ANDREW'S ALGORITHM RESULTANT CONVEX HULL :\n";
            for(i=0;i<Hull.size();i++)
                cout<<Hull[i].index<<" ("<<Hull[i].x<<", "<<Hull[i].y<<")"<<endl;
        }

        cout<<endl<<"Press any key to proceed."<<endl;
        getch();
        cout<<"\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
        cout<<"Enter 1-4 to choose which algorithm to use."<<endl;
        cout<<"1. Jarvis's March."<<endl;
        cout<<"2. Graham's Scan."<<endl;
        cout<<"3. Andrew's Algorithm."<<endl;
        cout<<"4. Use Visualizer for this result and Exit."<<endl;
        cout<<"5. Enter different points."<<endl;
        cout<<"6. Exit."<<endl;
        cin>>choice;

        if(choice==4){obj.toVisualizer(points,Hull,n);}
        if(choice==5){goto x;}

    }*/
	obj.toVisualizer(points,Hull,n);
    return 0;
}
