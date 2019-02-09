#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <functional>

//#####################################################
// FOLLOWING ARE DIFFERENT TEST CASES I DESIGNED
// SEE main() FOR DETAILS
// UNCOMMENT ANY DEFINE FOR TEST
//#####################################################

//Detailed DEBUG for collision
//#define COLLISION_DEBUG
//#define ROTATION_DEBUG -- For two scenarios, testing rotation of polygons


//DEBUG -- Turns on Warning messages about failed operations
//#define DEBUG

#define BOARD_DEBUG
#define PLAYER_DEBUG
#define PLAYER_DRAGGED_OUT_OF_BOARD
#define PLAYER_DRAGGED_COLLISION
#define PLAYER_ROTATED_OUT_OF_BOARD
#define PLAYER_ROTATED_COLLISION


using namespace std::placeholders;

//#####################################################
// UTILITIES
//#####################################################

class Shape;
using mType = double;

//Struct used to report the max and min coordinates
//of shapes in each direction -- Assumed 2D
struct MaxMinCoords{
    mType xMin, xMax, yMin, yMax;
};

//Data Type for points
struct Point2D{
    Point2D() : x(0), y(0) {}
    Point2D(mType xc, mType yc) : x(xc), y(yc) {}
    mType x, y;
};

//Dot product in 2D
inline mType vectorDotProd2D(const Point2D &p1,
                             const Point2D &p2){
    return p1.x*p2.x+p1.y*p2.y;
}

//Cross product in 2D
inline mType vectorXProd2D(const Point2D &p1,
                           const Point2D &p2){
    return p1.x*p2.y-p1.y*p2.x;
}

using mPointType = Point2D;
//using mPointType = Point3D; //Can be another option

class Matrix2D{
public:
    Matrix2D(const std::vector<mType> &inputMat) : matrix(inputMat) {};
    //Multiptly matrix by vector
    //return vector
    mPointType dot(const mPointType &p){
            return mPointType(matrix[0]*p.x+matrix[1]*p.y,
                              matrix[2]*p.x+matrix[3]*p.y);
    }
private:
    std::vector<mType> matrix;
};

//Operators are overloaded for the point type
mPointType operator+(const mPointType &p1, const mPointType &p2){
    return mPointType(p1.x+p2.x,p1.y+p2.y);
}

mPointType& operator+=(mPointType &p1, const mPointType &p2){
    p1.x += p2.x;
    p1.y += p2.y;
    return p1;
}

mPointType operator-(const mPointType &p1, const mPointType &p2){
    return mPointType(p1.x-p2.x,p1.y-p2.y);
}

mPointType& operator-=(mPointType &p1, const mPointType &p2){
    p1.x -= p2.x;
    p1.y -= p2.y;
    return p1;
}

std::ostream &operator<<(std::ostream &os, const mPointType &pnt){
    os << pnt.x << ' ' << pnt.y;
    return os;
}

//#####################################################
// SHAPE BASE CLASS AND DERIVED CLASSES
//#####################################################

//Forward Declaration of ShapeVisitor Class
class ShapeVisitor;

//Interface -- Base class of all available shapes
class Shape{
public:
    Shape(const mPointType &cPointCoords) : centerPoint(cPointCoords) {};
    virtual ~Shape(){};

    virtual void accept(ShapeVisitor &visitor) = 0;

    /*
     * Return the maximum and minimum points of
     * the smallest bounding rectangle
     * that the shape can fit into
     */
    virtual MaxMinCoords boundingBoxLimits() const = 0;

    //Check if the shape contains a given point
    virtual bool isInside(const mPointType &pnt) const = 0;

    mPointType centerPoint;
};

/*
 * Abstract class for shapes defined with Vertices in 2D
 * Number of vertices is required to construct
 */
class Polygon2D : public Shape{
public:
    Polygon2D(const mPointType &cPointCoords, const std::size_t &vertCnt) :
    Shape(cPointCoords), vertices(vertCnt,{0,0}) {};

    std::vector<mPointType> vertices;

    virtual ~Polygon2D(){};
    void accept(ShapeVisitor &visitor);

    MaxMinCoords boundingBoxLimits() const{
        MaxMinCoords mPnts;
        auto maxminX = std::minmax_element(vertices.begin(),vertices.end(),
                                           [](mPointType p1, mPointType p2)
                                           {return p1.x < p2.x;});
        auto maxminY = std::minmax_element(vertices.begin(),vertices.end(),
                                           [](mPointType p1, mPointType p2)
                                           {return p1.y < p2.y;});
        mPnts.xMin = maxminX.first->x;
        mPnts.xMax = maxminX.second->x;
        mPnts.yMin = maxminY.first->y;
        mPnts.yMax = maxminY.second->y;
        return mPnts;
    }
};

/*
 * Circle class needs radius as an input
 */
class Circle : public Shape{
public:
    Circle(const mPointType &cPointCoords,
           const mType &inputRadius) : Shape(cPointCoords),
           centerPoint(cPointCoords), radius(inputRadius) {}

    void accept(ShapeVisitor &visitor);

    MaxMinCoords boundingBoxLimits() const{
        MaxMinCoords mPnts;
        mPointType maxPnt = centerPoint+mPointType(radius,radius);
        mPointType minPnt = centerPoint-mPointType(radius,radius);
        mPnts.xMax = maxPnt.x;
        mPnts.yMax = maxPnt.y;
        mPnts.xMin = minPnt.x;
        mPnts.yMin = minPnt.y;
        return mPnts;
    }

    //Boundary is considered as inside
    bool isInside(const mPointType &pnt) const{
        mPointType relPnt = pnt - centerPoint;
        return relPnt.x*relPnt.x + relPnt.y*relPnt.y <= radius*radius;
    }
    mPointType centerPoint;
private:
    mType radius;
};

/*
 * Rectangle is a Polygon2D with 4 vertices
 * needs cPoint, width, height and rotation angle
 * rotation angle is to be provided with respect
 * to the center of the Rectangle
 */
class Rectangle : public Polygon2D{
public:
    Rectangle(const mPointType &cPointCoords,
              const mType &inputWidth, const mType &inputHeight,
              const mType &inputAngle) :
              Polygon2D(cPointCoords, 4), width(inputWidth),
			  height(inputHeight), rotAngle(inputAngle) {

        //Conversion from degrees to radian
        mType radRotAngle = rotAngle * M_PI/180.;
        //Rotation Matrix in 2D
        Matrix2D rot2DMat({cos(radRotAngle),sin(radRotAngle),
                          -sin(radRotAngle),cos(radRotAngle)});

        vertices[0] = rot2DMat.dot(cPointCoords +
                                   mPointType(-width/2.,-height/2.));
        vertices[1] = rot2DMat.dot(cPointCoords +
                                   mPointType(-width/2., height/2.));
        vertices[2] = rot2DMat.dot(cPointCoords +
                                   mPointType( width/2., height/2.));
        vertices[3] = rot2DMat.dot(cPointCoords +
                                   mPointType( width/2.,-height/2.));
    }

    /*
     * If a point(P) is inside the rectangle with edges (AB) and (AC)
     * then 0 < P.(AB) < |AB| and 0 < P.(AC) < |AC|
     * The edges are considered as inside
     */
    bool isInside(const mPointType &pnt) const{
        mPointType relPnt = pnt - vertices[0];
        mPointType edge1  = vertices[1]-vertices[0];
        mPointType edge2  = vertices[3]-vertices[0];

        mType projOnEdge1 = vectorDotProd2D(relPnt,edge1);
        mType edge1LenSqr = vectorDotProd2D(edge1,edge1);

        if(projOnEdge1 < 0 || projOnEdge1 > edge1LenSqr)
            return false;

        mType projOnEdge2 = vectorDotProd2D(relPnt,edge2);
        mType edge2LenSqr = vectorDotProd2D(edge2,edge2);

        if(projOnEdge2 < 0 || projOnEdge2 > edge2LenSqr)
            return false;

        return true;
    }
private:
    mType width, height, rotAngle;
};

/*
 * Triangle is a Polygon2D with 3 vertices
 * needs all the three vertices
 *
 * Sets the center point to -1,-1 in the initialization list
 * and then calculated within the constructor.
 *
 * ASSUMPTION: Input Triangle is a Triangle
 * I have not introduced any check mechanism for that
 */
class Triangle : public Polygon2D{
public:
    Triangle(const mPointType &vert1, const mPointType &vert2,
             const mPointType &vert3) :
    Polygon2D({-1,-1},3) {
        vertices[0] = vert1;
        vertices[1] = vert2;
        vertices[2] = vert3;
		mType cPointX = (vert2.x+vert3.x)/2.;
		cPointX += (vert1.x-cPointX)/3.;

		mType cPointY = (vert2.y+vert3.y)/2.;
		cPointY += (vert1.y-cPointY)/3.;

		//Set the centerPoint
		centerPoint = {cPointX,cPointY};
    }

    /*
     * Using three BaryCentric coordinates,
     * if the point is inside all must be >= 0
     * and the summation must be 1
     * (the Edges are considered inside)
     */
    bool isInside(const mPointType &pnt) const{
        /*
         * Area sign changes with the ordering of vertices
         * ClockWise vs Counter-ClockWise
         */
        mType Area = 0.5*(vectorXProd2D(vertices[1],vertices[2]) +
                          vectorXProd2D(vertices[2]-vertices[1],vertices[0]));

        mType s  = (vectorXProd2D(vertices[2],vertices[0]) +
                    vectorXProd2D(pnt,vertices[2]-vertices[0]))/(Area*2.);

        mType t  = (vectorXProd2D(vertices[0],vertices[1]) +
                    vectorXProd2D(pnt,vertices[0]-vertices[1]))/(Area*2.);

        return (s >= 0 && t >= 0 && 1-s-t >= 0);
    };
};

//#####################################################
//ShapeVisitor BASE CLASS AND DERIVED CLASSES
//#####################################################

//Interface -- Base class of all available visitors on Shapes
class ShapeVisitor{
public:
    virtual void execute(Circle    &inputShape) = 0;
    virtual void execute(Polygon2D &inputShape) = 0;
    virtual ~ShapeVisitor() {};
};

/*
 * DragVisitor
 * Input: move vector (Point2D type)
 * execute() method is overloaded for different shapes
 */
class DragVisitor : public ShapeVisitor{
public:
    DragVisitor(mPointType move) : moveVector(move) {};
    void execute(Circle &circle){
        //Update the center point of Circle
        circle.centerPoint += moveVector;
    }
    void execute(Polygon2D &poly){
        //Update the vertices
        std::for_each(poly.vertices.begin(),poly.vertices.end(),
                      [this](mPointType &p){p+=this->moveVector;});
    }
private:
    mPointType moveVector;
};
/*
 * RotateVisitor
 * Input: reference point           (Point2D type)
 *        rotation angle in degrees (Positive in Clock-Wise direction
 */
class RotateVisitor : public ShapeVisitor{
public:
    RotateVisitor(mPointType inputRefPoint, mType inputRotAngle) :
    refPoint(inputRefPoint), rotAngle(inputRotAngle) {};

    void execute(Circle &circle){
        //Update the center point of Circle
        rotatePoint2D(circle.centerPoint,refPoint);
    }
    void execute(Polygon2D &poly){
        //Update the vertices
    	/*
    	 * After rotating center point,
    	 * rotate the vertices
    	 */
    	mPointType oldCenter = poly.centerPoint;

    	/*
    	 * Rotate in three steps
    	 * 1- Update center point
    	 * 1- Move vertices
    	 * 2- Rotate vertices with respect to
    	 * the new centerPoint
    	 */
    	rotatePoint2D(poly.centerPoint,refPoint);
        std::for_each(poly.vertices.begin(),poly.vertices.end(),
                      [this,&poly,&oldCenter](mPointType& p)
					  {p += (poly.centerPoint-oldCenter);
					   rotatePoint2D(p,poly.centerPoint);});
    }
private:
    /*
     * Using Rotation Matrix in 2D rotate the point with respect to
     * a reference point.
     * newCoord = RefPoint + Rot_Matrix . (Pnt-RefPoint)
     */
    void rotatePoint2D(mPointType &pnt, const mPointType &refPnt){
        //Convert from degrees to radian
        mType radRotAngle = rotAngle * M_PI/180.;
        mPointType newCarCoord = pnt - refPnt;

        Matrix2D rot2DMat({cos(radRotAngle),sin(radRotAngle),
                          -sin(radRotAngle),cos(radRotAngle)});

        pnt = refPnt + rot2DMat.dot(newCarCoord);
    }

    mPointType refPoint;
    mType rotAngle;
};

//Define Shape::accept
//The best solution would be to move the declaration to
//a header file, however I am not doing it in this case
//to comply with the request of keeping all the code
//in a single file.
void Circle::accept(ShapeVisitor &visitor) {
    visitor.execute(*this);
}
void Polygon2D::accept(ShapeVisitor &visitor) {
    visitor.execute(*this);
}

//#############################################################
//FUNCTION CLASS CheckOverLap
//Checks for Overlapping Objects in two steps
//1- See if the bounding rectangles overlap (checkBoundingBoxes)
//2- If they do, then resample the combined bounding
//box areas of the two shapes with a high sampling
//rate in accordance with the acceptable error margin.
//Query if both the shapes share a commong point. checkDenseSampling
//#############################################################

/*
 * Function class that checks if there is an overlap between shapes
 * defined as binary_function to compare two shapes at a time
 */
class CheckOverLap : public std::binary_function<Shape *, Shape *, bool>{
public:
    bool operator()(const Shape *s1, const Shape *s2){
        if(checkBoundingBoxes(s1,s2)){
#ifdef COLLISION_DEBUG
            std::cout << "Possible overlap "
                     	 "Needs further investigation"
                     	 "This is done by checkDenseSampling routine\n";
#endif
            //Possible overlap -- Needs further investigation
            return checkDenseSampling(s1,s2);
        }
        return false;
    }
private:
    /*
     * Sampling distance
     * Should be compatible with the GameBoard dimensions
     */
    const mType pixSize = 0.1;

    /*
     *Check if the bounding boxes overlap
     */
    bool checkBoundingBoxes(const Shape *s1, const Shape *s2) const{
        MaxMinCoords mm1 = s1->boundingBoxLimits();
        MaxMinCoords mm2 = s2->boundingBoxLimits();

        mType d1x = mm2.xMin - mm1.xMax;
        mType d1y = mm2.yMin - mm1.yMax;
        mType d2x = mm1.xMin - mm2.xMax;
        mType d2y = mm1.yMin - mm2.yMax;

        if (d1x > 0.0 || d1y > 0.0)
			return false;

        if (d2x > 0.0 || d2y > 0.0)
			return false;

        return true;
    }

    /*
     * Check if the two shapes overlap by querying if points are
     * inside both of them.
     * The points are equally spaced, located at a distance defined by
     * variable pixSize, inside the bounding box covering the
     * two shapes.
     * MAX ERROR ~ pixSize
     */
    bool checkDenseSampling(const Shape *s1, const Shape *s2) const{
        MaxMinCoords mm1 = s1->boundingBoxLimits();
        MaxMinCoords mm2 = s2->boundingBoxLimits();

        //Combined Bounding Box
        mType xMin = std::min(mm1.xMin,mm2.xMin);
        mType xMax = std::max(mm1.xMax,mm2.xMax);

        mType yMin = std::min(mm1.yMin,mm2.yMin);
        mType yMax = std::max(mm1.yMax,mm2.yMax);


#ifdef COLLISION_DEBUG
        std::cout << xMin << ' ' << xMax << std::endl;
        std::cout << yMin << ' ' << yMax << std::endl;
#endif

        mType dX = xMax-xMin;
        mType dY = yMax-yMin;

        //Sampling points grid size
        unsigned nx = (unsigned)std::ceil(dX/pixSize);
        unsigned ny = (unsigned)std::ceil(dY/pixSize);

#ifdef COLLISION_DEBUG
        std::cout << dX << ' ' << dY << std::endl;
        std::cout << nx << ' ' << ny << std::endl;
#endif

        /*
         * QUERY if there is any cell shared by two shapes
         * Center and corner points of the imaginary grid cell
         * THE UGLIEST PART OF THE CODE IN MY OPINION
         */
        for(unsigned i = 0; i < nx; i++){
            for(unsigned j = 0; j < ny; j++){

                mPointType pnt1 = {xMin+(i+0.5)*pixSize,yMin+(j+0.5)*pixSize};
                mPointType pnt2 = {xMin+(i)*pixSize,yMin+(j)*pixSize};
                mPointType pnt3 = {xMin+(i)*pixSize,yMin+(j+1)*pixSize};
                mPointType pnt4 = {xMin+(i+1)*pixSize,yMin+(j)*pixSize};
                mPointType pnt5 = {xMin+(i+1)*pixSize,yMin+(j+1)*pixSize};

                if((s1->isInside(pnt1) && s2->isInside(pnt1)) ||
                   (s1->isInside(pnt2) && s2->isInside(pnt2)) ||
                   (s1->isInside(pnt3) && s2->isInside(pnt3)) ||
                   (s1->isInside(pnt4) && s2->isInside(pnt4)) ||
                   (s1->isInside(pnt5) && s2->isInside(pnt5))){
#ifdef COLLISION_DEBUG
                    std::cout << "Found shared point\n" << std::endl;
#endif
                    return true;
                }
            }
        }
        return false;
    }
};

//#####################################################
//GAME BOARD CLASS
//DRIVER OF THE GAME
//#####################################################

//A suggestion for GameBoard
class GameBoard{
public:
    GameBoard(mPointType inputBoardDims) :
    player(nullptr), boardDims(inputBoardDims){};

    ~GameBoard(){
        /*
         * Ensure that the shape pointers are deleted once the GameBoard
         * is out of scope.
         * This may not be necessary if the shapes are deleted by the
         * caller of GameBoard.
         * In this case, I am assuming that the shapes are provided as
         * pointers and the deletion is done by GameBoard.
         * See main() Test Simulation, to see how I am expecting it to
         * work.
         * I am not overriding copy and assignment operators not to dig
         * more into the implementation of GameBoard.
         * Otherwise, with destructor, copy and assignment operators should
         * be redefined as well.
         */
        std::for_each(boardShapes.begin(),boardShapes.end(),
                              [](Shape *s){delete s;});
        if(player != nullptr)
                delete player;
    }
    bool insertShape(Shape *newShape){
        if(!isInsideBoard(newShape)){
                std::cerr << "Shape is located out of the Board!";
                return false;
        }
        if(checkCollision(newShape)){
#ifdef DEBUG
            std::cerr << "New Shape cannot be Inserted! " <<
                         "-- Overlap with Board Shape\n";
#endif
            return false;
        }
        boardShapes.push_back(newShape);
        return true;
    }

    //Assuming memory is released after this operation
    bool removeShape(Shape *s){
        auto newEnd = std::remove(boardShapes.begin(),boardShapes.end(),s);
        if(newEnd != boardShapes.end()){
            boardShapes.erase(newEnd);
            return true;
        } else {
            std::cerr << "Shape is not on the board!\n";
            return false;
        }
    }

    bool setPlayer(Shape *myPlayer){
        if(player != nullptr){
#ifdef DEBUG
            std::cerr << "Already has a player\n";
#endif
            return false;
        }
        if(!isInsideBoard(myPlayer)){
#ifdef DEBUG
            std::cerr << "Player is located out of the Board!\n";
#endif
            return false;
        }
        if(checkCollision(myPlayer)){
#ifdef DEBUG
            std::cerr << "Player collides with Board Shape!\n";
#endif
            return false;
        }
        player = myPlayer;
        return true;
    }

    bool removePlayer(){
        if(player != nullptr){
            delete player;
            player = nullptr;
            return true;
        } else {
#ifdef DEBUG
            std::cerr << "Trying to remove non-existing Player!\n";
#endif
            return false;
        }
    }

    //Control the game through GameBoard by passing Commands
    bool passCommands(ShapeVisitor &visitor){
        player->accept(visitor);

        if(!isInsideBoard(player)){
#ifdef DEBUG
            std::cerr << "Player is out of the Board!\n";
#endif
            return false;
        }
        if(checkCollision(player)){
#ifdef DEBUG
            std::cerr << "Collision with Board Shape\n";
#endif
            return false;
        }
        return true;
    }
private:
    std::vector<Shape*> boardShapes;
    Shape *player;
    mPointType boardDims;

    bool isInsideBoard(const Shape *newShape) const{
        MaxMinCoords myBoundaries = newShape->boundingBoxLimits();
        if(myBoundaries.xMin < 0 || myBoundaries.xMax > boardDims.x ||
           myBoundaries.yMin < 0 || myBoundaries.yMax > boardDims.y){
        	return false;
        }
        return true;
    }

    bool checkCollision(const Shape *newShape) const{
	auto sIter = std::find_if(boardShapes.begin(),boardShapes.end(),
                              std::bind(CheckOverLap(),_1,newShape));
        if(sIter != boardShapes.end()){
            return true;
        }
        return false;
    }
};

//Test Simulation
/*
 * Assuming that the Simulator doesn't violate the following assumptions:
 * 1- It creates only one Player. I did not impose any restriction on it.
 * 2- It deletes the shapes, players that are rejected by the GameBoard.
 * 3- It does not delete the accepted players, shapes itself.
 *    GameBoard takes care of it.
 */

int main(){

    /* Create a GameBoard
     * Set Field Dimensions to 100 x 100 (Unitless)
     * GameBoard Dimensions are used to check if shapes are
     * placed or moved out of the board.
     */
    GameBoard myGameBoard({100,100});

    Shape *bCircle = new Circle({30,30},5);
    Shape *bRect   = new Rectangle({70,70},10,20,0);
    Shape *bTri    = new Triangle({50,50},{50,60},{60,50});

    if(!myGameBoard.insertShape(bCircle)){
        delete bCircle;
    }
    if(!myGameBoard.insertShape(bRect)){
         delete bRect;
    }
    if(!myGameBoard.insertShape(bTri)){
        delete bTri;
    }

#ifdef ROTATION_DEBUG

    {
		Triangle *myPlayer = new Triangle({0,0},{0,20},{10,0});
		RotateVisitor   rot({0,0},90);

		std::cout << "----- Rotating Triangle by 90 deg\n";
		std::cout << "Original Location: \n";
		std::cout << "\tCenter Point: " << myPlayer->centerPoint << std::endl;
		std::cout << "Vertices:\n";
		unsigned int i = 1;
		for(auto it : myPlayer->vertices){
			std::cout << "\tVertex - " << i++ << ": " << it << std::endl;
		}

		if(!myGameBoard.setPlayer(myPlayer)){
			delete myPlayer;
		}

		if(!myGameBoard.passCommands(rot)){
			std::cerr << "Could not pass the command!!!\n";
		}

		std::cout << std::endl;

		std::cout << "New Location: \n";
		std::cout << "\tCenter Point: " << myPlayer->centerPoint << std::endl;
		std::cout << "Vertices:\n";
		i = 1;
		for(auto it : myPlayer->vertices){
			std::cout << "\tVertex - " << i++ << ": " << it << std::endl;
		}

		myGameBoard.removePlayer();
    }

    std::cout << std::endl;

    {

		Rectangle *myPlayer = new Rectangle({10,10},10,10,0);
		RotateVisitor   rot({0,0},45);

		std::cout << "----- Rotating Rectangle by 45 deg\n";
		std::cout << "Original Location: \n";
		std::cout << "\tCenter Point: " << myPlayer->centerPoint << std::endl;
		std::cout << "Vertices:\n";
		unsigned int i = 1;
		for(auto it : myPlayer->vertices){
			std::cout << "\tVertex - " << i++ << ": " << it << std::endl;
		}

		if(!myGameBoard.setPlayer(myPlayer)){
			delete myPlayer;
		}

		if(!myGameBoard.passCommands(rot)){
			std::cerr << "Could not pass the command!!!\n";
		}

		std::cout << std::endl;

		std::cout << "New Location: \n";
		std::cout << "\tCenter Point: " << myPlayer->centerPoint << std::endl;
		std::cout << "Vertices:\n";
		i = 1;
		for(auto it : myPlayer->vertices){
			std::cout << "\tVertex - " << i++ << ": " << it << std::endl;
		}
    }
#endif

#ifdef BOARD_DEBUG
    {
        std::cout << "########BOARD DEBUG###############################\n";
        /*
         * Overlapping Circles that can be identified
         * by refined point search -- see checkDenseSampling()
         */
        //Circle({30,30},5) - Is already on the Board
        Shape *c1 = new Circle({22,22},5);
        if(myGameBoard.insertShape(c1)){ //Can be inserted
            std::cout << "TEST SUCCEED -- Inserted Board Shape\n";
            //Removing shape not to cause problem with other shapes
            myGameBoard.removeShape(c1);
            delete c1;
        } else {
            std::cout << "TEST FAILED"
                         " -- Cannot Insert Non-Overlapping Shape\n";
            delete c1;
        }
        Shape *c2 = new Circle({37,37},5); //Overlaps
        if(myGameBoard.insertShape(c2)){   //Cannot be inserted
            std::cout << "TEST FAILED -- Inserted Board Shape Overlapping\n";
        } else {
            std::cout << "TEST SUCCEED -- Cannot Insert Overlapping Shape"
                         " -- Circle over Circle\n";
            delete c2;
        }

        Shape *r1 = new Rectangle({37,37},14,14,0); //Overlaps
        if(myGameBoard.insertShape(r1)){   //Cannot be inserted
            std::cout << "TEST FAILED -- Inserted Board Shape Overlapping\n";
        } else {
            std::cout << "TEST SUCCEED -- Cannot Insert Overlapping Shape"
                         " -- Rectangle over Circle\n";
            delete r1;
        }

        Shape *t1 = new Triangle({32,32},{42,42},{36,50}); //Overlaps
        if(myGameBoard.insertShape(t1)){   //Cannot be inserted
            std::cout << "TEST FAILED -- Inserted Board Shape Overlapping\n";
        } else {
            std::cout << "TEST SUCCEED -- Cannot Insert Overlapping Shape"
                         " -- Triangle over Circle\n";
            delete t1;
        }

        //new Triangle({50,50},{50,60},{60,50}) already on board
        Shape *t2 = new Triangle({55,55},{60,60},{60,70}); //Overlaps
        if(myGameBoard.insertShape(t2)){   //Cannot be inserted
            std::cout << "TEST FAILED -- Inserted Board Shape Overlapping\n";
        } else {
            std::cout << "TEST SUCCEED -- Cannot Insert Overlapping Shape"
                         " -- Triangle over Triangle\n";
            delete t2;
        }

        std::cout << "########BOARD DEBUG###############################\n";
        std::cout << std::endl;
    }
#endif

#ifdef PLAYER_DEBUG
    {
        std::cout << "########PLAYER DEBUG##############################\n";

        //Player overlaps with Board Shape
        Shape *myPlayer1 = new Circle({10,10},5);
        if(myGameBoard.setPlayer(myPlayer1)){ //Success
            std::cout << "TEST SUCCEED -- Set Non-Overlapping Player\n";
        } else {
            std::cout << "TEST FAILED -- Set Overlapping Player\n";
            delete myPlayer1;
        }
        myGameBoard.removePlayer();

        Shape *myPlayer2 = new Triangle({30,30},{30,40},{40,30});
        if(myGameBoard.setPlayer(myPlayer2)){ //Failure
            std::cout << "TEST FAILED -- Set the Overlapping Player\n";
            myGameBoard.removePlayer();
        } else {
            std::cout << "TEST SUCCEED -- Cannot set the Overlapping Player\n";
            delete myPlayer2;
        }
        std::cout << "########PLAYER DEBUG##############################\n";
        std::cout << std::endl;
    }
#endif

#ifdef PLAYER_DRAGGED_OUT_OF_BOARD
    {
        std::cout << "########PLAYER DRAG OUT OF BOARD DEBUG############\n";

        Shape *myPlayer = new Circle({10,10},5);
        DragVisitor   drag1({1,0});
        DragVisitor   drag2({-11,0});

        if(!myGameBoard.setPlayer(myPlayer)){
            delete myPlayer;
        }
        if(myGameBoard.passCommands(drag1)){ //Success
            std::cout << "TEST SUCCEED -- Dragged the Player on Board\n";
        } else {
            std::cout << "TEST FAILED -- Warning about Dragging "
                         "the Player on Board -- No warning expected\n";
        }

        if(myGameBoard.passCommands(drag2)){ //Failure
            std::cout << "TEST FAILED -- No Warning about"
                         "Dragged the Player out of Board"
                         " -- No warning expected\n";
        } else {
            std::cout << "TEST SUCCEED -- Warning about "
                         "Dragging the Player out of Board\n";
        }

        myGameBoard.removePlayer();

        std::cout << "########PLAYER DRAG OUT OF BOARD DEBUG############\n";
        std::cout << std::endl;
    }
#endif

#ifdef PLAYER_DRAGGED_COLLISION
    {
        std::cout << "########PLAYER DRAGGED COLLISION DEBUG############\n";

        Shape *myPlayer = new Rectangle({10,10},20,20,0);
        DragVisitor   drag1({1,1});
        DragVisitor   drag2({49,49});

        if(!myGameBoard.setPlayer(myPlayer)){
            delete myPlayer;
        }

        if(myGameBoard.passCommands(drag1)){ //Success
             std::cout << "TEST SUCCEED -- Dragged the Player"
                          " -- No Collision Expected\n";
        } else {
            std::cout << "TEST FAILED -- Warning about Dragging the Player"
                         " -- No Collision Expected\n";
        }
        if(myGameBoard.passCommands(drag2)){ //Failure
            std::cout << "TEST FAILED -- Dragged the Player to Collision"
                         " -- Collision Expected\n";
        } else {
            std::cout << "TEST SUCCEED -- Warning about Dragging the Player"
                         " -- Collision Expected\n";
        }

        myGameBoard.removePlayer();
        std::cout << "########PLAYER DRAGGED COLLISION DEBUG############\n";
        std::cout << std::endl;
    }
#endif

#ifdef PLAYER_ROTATED_OUT_OF_BOARD
    {
        std::cout << "########PLAYER ROTATE OUT OF BOARD DEBUG##########\n";
        /*
         * 14.04 deg CW rotation with respect to board {0,0}
         * moves it out of board.
         */
        Shape *myPlayer = new Rectangle({10,15},20,20,0);
        RotateVisitor   rot1({0,0},14);  //--SAFE
        RotateVisitor   rot2({0,0},-15); //--FAIL

        if(!myGameBoard.setPlayer(myPlayer)){
            delete myPlayer;
        }
        if(myGameBoard.passCommands(rot1)){
            std::cout << "TEST SUCCEED -- Rotated the Player\n";
        } else{
            std::cout << "TEST FAILED -- Warning about Rotating "
                         "the Player -- No Warning Expected\n";
        }

        if(myGameBoard.passCommands(rot2)){
            std::cout << "TEST FAILED -- No Warning about "
                         "moving Out of Board\n";
        } else {
            std::cout << "TEST SUCCEED -- Warning about moving"
                         " Out of Board\n";
        }

        myGameBoard.removePlayer();
        std::cout << "########PLAYER ROTATE OUT OF BOARD DEBUG##########\n";
        std::cout << std::endl;
    }
#endif

#ifdef PLAYER_ROTATED_COLLISION
    {
        std::cout << "########PLAYER ROTATE TO COLLISION DEBUG##########\n";
        Circle *myPlayer = new Circle({35,24},2);
        RotateVisitor   rot1({0,0},1);   //--SAFE
        RotateVisitor   rot2({0,0},-20); //--FAIL

        if(!myGameBoard.setPlayer(myPlayer)){
            delete myPlayer;
        }
        if(myGameBoard.passCommands(rot1)){ //Success
            std::cout << "TEST SUCCEED -- Rotated the Player\n";
        } else {
            std::cout << "TEST FAILED -- Cannot Rotate the Player\n";
        }
        if(myGameBoard.passCommands(rot2)){ //Failure
            std::cout << "TEST FAILED -- No Warning about Collision"
                         "Rotated the Player\n";
        } else {
            std::cout << "TEST SUCCEED -- Warning about Collision\n";
        }

        myGameBoard.removePlayer();
        std::cout << "########PLAYER ROTATE TO COLLISION DEBUG##########\n";
        std::cout << std::endl;
    }
#endif
}

