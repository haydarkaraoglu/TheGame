# The Game

## Description of Task

You are part of a team that develops a small children's computer game. In this game a number of objects are shown on the screen. Each object consists of a number of simple geometric shapes: circles, triangles and rectangles. 

One of the objects is the player's playing piece. His task is to move his playing piece from its starting position to a set end position without bumping into any of the other objects. He can drag his playing piece or rotate it about a chosen point.

__Your task is to model the game objects and implement the math. Write functions that implement the operations that the player can do and to detect collisions.__

## Overview


## Build steps

Program is written in C++11. As compiler I used *g++ (Ubuntu 7.3.0-27ubuntu1~18.04) 7.3.0*. To build run:

```
g++ theGame.cpp
```

## Sample execution

```bash
./theGame #ifdef PLAYER_DEBUG
TEST SUCCEED -- Cannot set the Overlapping Player
TEST SUCCEED -- Set Non-Overlapping Player
...

```

## Architecure


#### PATTERN:
* I followed VISITOR DESIGN PATTERN to separate operations from shapes.
* Shape accepts ShapeVisitor passing its reference to it and ShapeVisitor updates Shape by executing its tasks.

#### SHAPE:
* Both Rectangle and Triangle are Polygons in 2D. Thus, I have an abstract class Polygon2D that defines a member function (boundingBoxLimits()) which is the same for both Rectangle and Triangle.

	___IDEA___: To add other curvy shapes such as ellipse, it might be an option to add another Abstract Class to refrain from code repetition.

* I made the centerPoint for Circles and vertices for Polygon2Ds public, so that visitors can access and modify them. 

#### OPERATION:

1- Classes using ShapeVisitor Interface execute their tasks on the Shapes they receive. I have two execute() routines, one for Circle and another for Polygon2D type objects. This class needs to be extended for every new shape to be introduced.

2- DragVisitor moves the Shape.

3- RotateVisitor rotates Shape with respect to a reference point and for a given angle (positive in Clock-wise direction - in degrees).

#### CHECKOVERLAP:

1. CheckOverLap functor (inheriting from binary_function class) compares two shapes locations in two steps. 

	1. First, minimum bounding boxes are compared based on the bounding box limits returned from shapes. This is a quick check and eliminates shape pairs that are located far enough from each other. 
	2. If two bounding boxes overlap, the functor calls another member function that samples a new minimum bounding box covering the two shapes. The sampling rate is defined by a private attribute (pixSize). The error in the detection is around pixSize.
 
	___GENERAL COMMENTS ON CHECKOVERLAP___:

	Another option to check collision would be through geometrical analysis. This would be exact, however such an approach would require modifications of CheckOverLap class each time a new shape is added for different combinations of shapes. Also for complex shapes such as concave ones, I expect geometrical analysis to be quite demanding. 

	Thus, I preferred to come up with a more general solution that relies on sampling the minimum bounding box. Choosing a sampling rate in accordance with the expected resolution is preferrable. This approach requires each shape to provide a method to check if a point is inside or outside of it, thus delegating the responsibility to shape objects. I implemented isInside() method for the three shapes. Details and the underlying math are explained in the code as comments.

#### GAMEBOARD:

To my understanding, this is not really a part of the assignment, however I implemented a GameBoard for testing purposes:

1. GameBoard is the controller of the game. It has dimensions and the main coordinate system is defined on it with origin assumed to be located at (0,0) in 2D.
2. It sets a single player (shape) and keeps record of the other shapes located on the board.
3. For each operation (visitor), it checks if the player shape stays in the Board and if there is any collision.

#### UTILITIES:

1. I defined a struct to use as my point data type.
2. Some operators and vector algebra tools are provided to facilitate the implementation.

#### ASSUMPTIONS:
1. Triangle: I assume the provided vertices form a triangle.
2. The shapes inserted in the GameBoard once accepted belong to GameBoard and it deletes them upon completion. If the shape is rejected it is the responsibility of the user of the GameBoard to take care of deletion. 
 
	__IDEA__: Using smart pointers would be a better option here, but as designing GameBoard is out of the scope of this task, I did not spend much time on it.

