#include "oglwidget.h"
#include <math.h>
#include <iostream>
#include <QDebug>
#include <complex>
#include <deque>

#define PI 3.14159265358979323846
using namespace std;

static double alpha = 45.0; // rotation angle
static double beta = 30.0; // rotation angle
static double gamma3 = 85.0; // rotation angle
static double angleTorus1 = 0.0;
static double angleTorus2 = 60.0;
static double angleTorus3 = 120.0;
static double angleTorus4 = 180.0;
static double angleTorus5 = 240.0;
static double angleTorus6 = 300.0;

static float doublePI = 2 * PI;


// initialize Open GL lighting and projection matrix
void InitLightingAndProjection() // to be executed once before drawing
{
    // light positions and colors
    GLfloat LightPosition1[4] = { 10, 5, 10,  0};
    GLfloat LightPosition2[4] = { -5, 5, -10,  0};
    GLfloat ColorRedish[4] = { 1.0,  .8,  .8,  1}; // white with a little bit of red
    GLfloat ColorBlueish[4] = { .8,  .8,  1.0,  1};// white with a little bit of blue

    glEnable( GL_DEPTH_TEST); // switch on z-buffer
    glDepthFunc( GL_LESS);

    glShadeModel( GL_SMOOTH); // Gouraud shading
    //glShadeModel( GL_FLAT);

    glEnable( GL_LIGHTING); // use lighting
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 1); // draw both sides

    // define and switch on light 0
    glLightfv( GL_LIGHT0, GL_POSITION, LightPosition1);
    glLightfv( GL_LIGHT0, GL_DIFFUSE,  ColorRedish);
    glLightfv( GL_LIGHT0, GL_SPECULAR, ColorRedish);
    glEnable( GL_LIGHT0);

    // define and switch on light 1
    glLightfv( GL_LIGHT1, GL_POSITION, LightPosition2);
    glLightfv( GL_LIGHT1, GL_DIFFUSE,  ColorBlueish);
    glLightfv( GL_LIGHT1, GL_SPECULAR, ColorBlueish);
    glEnable( GL_LIGHT1);

    glMatrixMode( GL_PROJECTION); // define camera projection
    glLoadIdentity(); // reset matrix to identity (otherwise existing matrix will be multiplied with)
    glOrtho( -15, 15, -10, 10, -20, 20); // orthogonal projection (xmin xmax ymin ymax zmin zmax)
    //glFrustum( -10, 10, -8, 8, 2, 20); // perspective projektion
}

void DrawCylinder( int reso = 16){ // drawing a cylinder in OpenGL
    // alocate memory for x and y coordinates on a circle
    double *c = new double[ reso+1];
    double *s = new double[ reso+1];

    for( int i=0; i<=reso; i++){ // compute x and y coordinates of citcle
        c[i] = cos( 2.0 * PI * i / reso );
        s[i] = sin( 2.0 * PI * i / reso );
        //cout << i << " " << c[i] << endl;
    }

    glBegin( GL_QUADS); // each 4 points define a polygon
    for( int i=0; i<reso; i++){
        glNormal3f( c[i], s[i], 0.0); // normal vector used for all consecutive points
        glVertex3f( c[i], s[i], 3.0); // 2 points ...
        glVertex3f( c[i], s[i], 0.0);

        glNormal3f( c[i+1], s[i+1], 0.0); // another normal with two more points
        glVertex3f( c[i+1], s[i+1], 0.0);
        glVertex3f( c[i+1], s[i+1], 3.0);
    }
    glEnd(); // concludes GL_QUADS

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 0);
    glVertex3f(0, 0, 0); // center of circle
    for( int j = 0; j <= reso;j++) {
        glVertex3f(
                    sin( 2.0 * PI * j / reso ),
                    cos( 2.0 * PI * j / reso ),
                    0
                    );

    }
    glEnd();
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 3);
    glVertex3f(0, 0, 3); // center of circle
    for( int j = 0; j <= reso;j++) {
        glVertex3f(
                    cos( 2.0 * PI * j / reso ),
                    sin( 2.0 * PI * j / reso ),
                    3
                    );

    }
    glEnd();

    delete[] c; // de-allocate space
    delete[] s;
}


void DrawCube(){ // drawing a cylinder in OpenGL

    int lenghtOfCube = 1;

    glBegin( GL_QUADS);

    //Z ist die nach vorne gerichtetet Seite.
    glNormal3f(0,1,0);//oben
    glVertex3f( lenghtOfCube, lenghtOfCube,-lenghtOfCube);    // Top Right Of The Quad (Top)
    glVertex3f(-lenghtOfCube, lenghtOfCube,-lenghtOfCube);    // Top Left Of The Quad (Top)
    glVertex3f(-lenghtOfCube, lenghtOfCube, lenghtOfCube);    // Bottom Left Of The Quad (Top)
    glVertex3f( lenghtOfCube, lenghtOfCube, lenghtOfCube);    // Bottom Right Of The Quad (Top)

    glNormal3f(0,-1,0);//unten
    glVertex3f( lenghtOfCube, -lenghtOfCube, lenghtOfCube);    // Top Right Of The Quad (Bottom)
    glVertex3f(-lenghtOfCube, -lenghtOfCube, lenghtOfCube);    // Top Left Of The Quad (Bottom)
    glVertex3f(-lenghtOfCube, -lenghtOfCube,-lenghtOfCube);    // Bottom Left Of The Quad (Bottom)
    glVertex3f( lenghtOfCube, -lenghtOfCube,-lenghtOfCube);    // Bottom Right Of The Quad (Bottom)

    glNormal3f(0,0,1);//vorne
    glVertex3f( lenghtOfCube, lenghtOfCube, lenghtOfCube);    // Top Right Of The Quad (Front)
    glVertex3f(-lenghtOfCube, lenghtOfCube, lenghtOfCube);    // Top Left Of The Quad (Front)
    glVertex3f(-lenghtOfCube,-lenghtOfCube, lenghtOfCube);    // Bottom Left Of The Quad (Front)
    glVertex3f( lenghtOfCube,-lenghtOfCube, lenghtOfCube);    // Bottom Right Of The Quad (Front)

    glNormal3f(0,0,-1);//hinten
    glVertex3f( lenghtOfCube,-lenghtOfCube, -lenghtOfCube);    // Top Right Of The Quad (Back)
    glVertex3f(-lenghtOfCube,-lenghtOfCube, -lenghtOfCube);    // Top Left Of The Quad (Back)
    glVertex3f(-lenghtOfCube, lenghtOfCube, -lenghtOfCube);    // Bottom Left Of The Quad (Back)
    glVertex3f( lenghtOfCube, lenghtOfCube, -lenghtOfCube);    // Bottom Right Of The Quad (Back)

    glNormal3f(-1,0,0);//links
    glVertex3f(-lenghtOfCube, lenghtOfCube, lenghtOfCube);    // Top Right Of The Quad (Left)
    glVertex3f(-lenghtOfCube, lenghtOfCube,-lenghtOfCube);    // Top Left Of The Quad (Left)
    glVertex3f(-lenghtOfCube,-lenghtOfCube,-lenghtOfCube);    // Bottom Left Of The Quad (Left)
    glVertex3f(-lenghtOfCube,-lenghtOfCube, lenghtOfCube);    // Bottom Right Of The Quad (Left)

    glNormal3f(1,0,0);//rechts
    glVertex3f( lenghtOfCube, lenghtOfCube,-lenghtOfCube);    // Top Right Of The Quad (Right)
    glVertex3f( lenghtOfCube, lenghtOfCube, lenghtOfCube);    // Top Left Of The Quad (Right)
    glVertex3f( lenghtOfCube,-lenghtOfCube, lenghtOfCube);    // Bottom Left Of The Quad (Right)
    glVertex3f( lenghtOfCube,-lenghtOfCube,-lenghtOfCube);    // Bottom Right Of The Quad (Right) - See more at: http://www.codemiles.com/c-opengl-examples/draw-3d-cube-using-opengl-t9018.html#sthash.179MIp09.dpuf
    glEnd(); // concludes GL_QUADS
}

void DrawPyramid(){ // drawing a cylinder in OpenGL

    int lenghtOfQuad = 1;

    glBegin( GL_QUADS); // each 4 points define a polygon

    //Boden
    glNormal3f(0,0,0);
    glVertex3f( lenghtOfQuad,0, lenghtOfQuad);    // Top Right Of The Quad (Bottom)
    glVertex3f( lenghtOfQuad,0,-lenghtOfQuad);    // Bottom Right Of The Quad (Bottom)
    glVertex3f(-lenghtOfQuad,0,-lenghtOfQuad);    // Bottom Left Of The Quad (Bottom)
    glVertex3f(-lenghtOfQuad,0, lenghtOfQuad);    // Top Left Of The Quad (Bottom)
    glEnd(); // concludes GL_QUADS

    int lenghtOfTriangle = 1;


    glBegin( GL_TRIANGLES);

    //hinten
    glNormal3f(0,.5,-.5);
    glVertex3f(lenghtOfTriangle,0,-lenghtOfTriangle);
    glVertex3f(-lenghtOfTriangle,0,-lenghtOfTriangle);
    glVertex3f(0,lenghtOfTriangle,0);

    //rechts
    glNormal3f(.5,.5,0);
    glVertex3f(lenghtOfTriangle,0,lenghtOfTriangle);
    glVertex3f(lenghtOfTriangle,0,-lenghtOfTriangle);
    glVertex3f(0,lenghtOfTriangle,0);

    //vorne
    glNormal3f(0,.5,.5);
    glVertex3f(-lenghtOfTriangle,0,lenghtOfTriangle);
    glVertex3f(lenghtOfTriangle,0,lenghtOfTriangle);
    glVertex3f(0,lenghtOfTriangle,0);

    //links
    glNormal3f(-.5,.5,0);
    glVertex3f(-lenghtOfTriangle,0,-lenghtOfTriangle);
    glVertex3f(-lenghtOfTriangle,0,lenghtOfTriangle);
    glVertex3f(0,lenghtOfTriangle,0);
    glEnd();

}

void DrawTorus(float r, float R){
    int reso = 40;
    float *s = new float[ reso+1];
    //vector < vector <float> > x;
    //vector < vector <float> > y;
    //vector < vector <float> > z;

    float x[reso + 1][reso + 1];
    float y[reso + 1][reso + 1];
    float z[reso + 1][reso + 1];
    float n_x[reso + 1][reso + 1];
    float n_y[reso + 1][reso + 1];
    float n_z[reso + 1][reso + 1];

    //deque < float > norm_d;


    float sx = .0;
    float sy = .0;
    float sz = .0;
    float tx = .0;
    float ty = .0;
    float tz = .0;

    float n[3];


    for( int i=0; i <= reso; i++){ // compute x and y coordinates of citcle
        s[i] = doublePI / reso * i;
        //cout << i << " " << c[i] << endl;
    }

    for ( int i = 0; i <= reso; i++){
        for ( int j = 0; j <= reso; j++){
            x[i][j] = cosf(s[i]) * (R + r * cosf(s[j]));
            y[i][j] = sinf(s[i]) * (R + r * cosf(s[j]));
            z[i][j] = r * sinf(s[j]);

            sx = -sinf(s[i]) * (R + r * cosf(s[j]));
            sy = cosf(s[i]) * (R + r * cosf(s[j]));
            sz = 0;

            tx = cosf(s[i]) * r * -sinf(s[j]);
            ty = sinf(s[i]) * r * -sinf(s[j]);
            tz = r * cosf(s[j]);


            n[0] = sy * tz - sz * ty;
            n[1] = sz * tx - sx * tz;
            n[2] = sx * ty - sy * tx;
            /*

            float calc = .0;
            for (int m = 0; m < 3; m++){
                calc += n[m] * n[m];
            }

            float norm = sqrtf(calc);
*/
            n_x[i][j] = n[0];// norm;
            n_y[i][j] = n[1];// norm;
            n_z[i][j] = n[2];// norm;

            /*
            x[i][j] += n[0];
            y[i][j] += n[1];
            z[i][j] += n[2];
*/

        }
    }

    glBegin(GL_QUADS);

    //float norbert[3];

    //for (int m = 0; m < 3; m++){
    //    norbert[m] = n[m];
    //}

    for (int k = 0; k < reso; k++){
        for (int l = 0; l < reso; l++){
            glNormal3f(n_x[l][k],n_y[l][k],n_z[l][k]);
            glVertex3f(x[l][k],y[l][k],z[l][k]);
            glNormal3f(n_x[l+1][k],n_y[l+1][k],n_z[l+1][k]);
            glVertex3f(x[l+1][k],y[l+1][k],z[l+1][k]);
            glNormal3f(n_x[l+1][k+1],n_y[l+1][k+1],n_z[l+1][k+1]);
            glVertex3f(x[l+1][k+1],y[l+1][k+1],z[l+1][k+1]);
            glNormal3f(n_x[l][k+1],n_y[l][k+1],n_z[l][k+1]);
            glVertex3f(x[l][k+1],y[l][k+1],z[l][k+1]);

        }
    }

    glEnd(); // concludes GL_QUADS

    delete[] s; // de-allocate space
}

void DrawSphere(float r, float R){
    int reso = 360;
    float *s = new float[ reso+1];
    //vector < vector <float> > x;
    //vector < vector <float> > y;
    //vector < vector <float> > z;

    float x[reso + 1][reso + 1];
    float y[reso + 1][reso + 1];
    float z[reso + 1][reso + 1];
    float n_x[reso + 1][reso + 1];
    float n_y[reso + 1][reso + 1];
    float n_z[reso + 1][reso + 1];

    //deque < float > norm_d;


    float sx = .0;
    float sy = .0;
    float sz = .0;
    float tx = .0;
    float ty = .0;
    float tz = .0;

    float n[3];


    for( int i=0; i <= reso; i++){ // compute x and y coordinates of citcle
        s[i] = doublePI / reso * i;
        //cout << i << " " << c[i] << endl;
    }

    for ( int i = 0; i <= reso; i++){
        for ( int j = 0; j <= reso; j++){
            x[i][j] = cosf(s[i]) * (R + r * cosf(s[j]));
            y[i][j] = sinf(s[i]) * (R + r * cosf(s[j]));
            z[i][j] = r * sinf(s[j]);

            sx = -sinf(s[i]) * (R + r * cosf(s[j]));
            sy = cosf(s[i]) * (R + r * cosf(s[j]));
            sz = 0;

            tx = cosf(s[i]) * r * -sinf(s[j]);
            ty = sinf(s[i]) * r * -sinf(s[j]);
            tz = r * cosf(s[j]);


            n[0] = sy * tz - sz * ty;
            n[1] = sz * tx - sx * tz;
            n[2] = sx * ty - sy * tx;
            /*

            float calc = .0;
            for (int m = 0; m < 3; m++){
                calc += n[m] * n[m];
            }

            float norm = sqrtf(calc);
*/
            n_x[i][j] = n[0];// norm;
            n_y[i][j] = n[1];// norm;
            n_z[i][j] = n[2];// norm;



        }
    }

    glBegin(GL_QUADS);

    for (int k = 0; k < reso / 4; k++){
        for (int l = 0; l < reso / 2; l++){
            glNormal3f(n_x[l][k],n_y[l][k],n_z[l][k]);
            glVertex3f(x[l][k],y[l][k],z[l][k]);
            glNormal3f(n_x[l+1][k],n_y[l+1][k],n_z[l+1][k]);
            glVertex3f(x[l+1][k],y[l+1][k],z[l+1][k]);
            glNormal3f(n_x[l+1][k+1],n_y[l+1][k+1],n_z[l+1][k+1]);
            glVertex3f(x[l+1][k+1],y[l+1][k+1],z[l+1][k+1]);
            glNormal3f(n_x[l][k+1],n_y[l][k+1],n_z[l][k+1]);
            glVertex3f(x[l][k+1],y[l][k+1],z[l][k+1]);
        }
    }


    for (int k = 0; k < reso / 4; k++){
        for (int l = reso / 2; l < reso; l++){
            glNormal3f(n_x[l][k],n_y[l][k],n_z[l][k]);
            glVertex3f(x[l][k],y[l][k],z[l][k]);
            glNormal3f(n_x[l+1][k],n_y[l+1][k],n_z[l+1][k]);
            glVertex3f(x[l+1][k],y[l+1][k],z[l+1][k]);
            glNormal3f(n_x[l+1][k+1],n_y[l+1][k+1],n_z[l+1][k+1]);
            glVertex3f(x[l+1][k+1],y[l+1][k+1],z[l+1][k+1]);
            glNormal3f(n_x[l][k+1],n_y[l][k+1],n_z[l][k+1]);
            glVertex3f(x[l][k+1],y[l][k+1],z[l][k+1]);

        }
    }


    for (int k = reso / 4 * 3; k < reso ; k++){
        for (int l = 0; l < reso / 2; l++){
            glNormal3f(n_x[l][k],n_y[l][k],n_z[l][k]);
            glVertex3f(x[l][k],y[l][k],z[l][k]);
            glNormal3f(n_x[l+1][k],n_y[l+1][k],n_z[l+1][k]);
            glVertex3f(x[l+1][k],y[l+1][k],z[l+1][k]);
            glNormal3f(n_x[l+1][k+1],n_y[l+1][k+1],n_z[l+1][k+1]);
            glVertex3f(x[l+1][k+1],y[l+1][k+1],z[l+1][k+1]);
            glNormal3f(n_x[l][k+1],n_y[l][k+1],n_z[l][k+1]);
            glVertex3f(x[l][k+1],y[l][k+1],z[l][k+1]);

        }
    }


    for (int k = reso / 4 * 3; k < reso ; k++){
        for (int l = reso / 2; l < reso; l++){
            glNormal3f(n_x[l][k],n_y[l][k],n_z[l][k]);
            glVertex3f(x[l][k],y[l][k],z[l][k]);
            glNormal3f(n_x[l+1][k],n_y[l+1][k],n_z[l+1][k]);
            glVertex3f(x[l+1][k],y[l+1][k],z[l+1][k]);
            glNormal3f(n_x[l+1][k+1],n_y[l+1][k+1],n_z[l+1][k+1]);
            glVertex3f(x[l+1][k+1],y[l+1][k+1],z[l+1][k+1]);
            glNormal3f(n_x[l][k+1],n_y[l][k+1],n_z[l][k+1]);
            glVertex3f(x[l][k+1],y[l][k+1],z[l][k+1]);

        }
    }



    glEnd(); // concludes GL_QUADS

    delete[] s; // de-allocate space
}

void DrawMoebius(){
    int reso = 50;
    vector < float > a;
    vector < float > r;

    float x[reso + 1][reso + 1];
    float y[reso + 1][reso + 1];
    float z[reso + 1][reso + 1];

    for (int i = 0; i <= reso; i++){
        a.push_back(doublePI / reso * i);
    }

    for (float i = -1.; i <= 1.1; i += 2. / (float)reso){
        if ( i >= float(1)){
            r.push_back(float(1));
            i = 1.2;
        }
        else{
            r.push_back(i);
        }
    }

    for (int i = 0; i <= reso; i++){
        for ( int j = 0; j <= reso; j++){
            x[i][j] = cosf(a[i]) * (1 + (r[j]/2) * cosf(a[i] / 2));
            y[i][j] = sinf(a[i]) * (1 + (r[j]/2) * cosf(a[i] / 2));
            z[i][j] = (r[j] / 2) * sinf(a[i] / 2);
        }
    }

    glBegin(GL_QUADS);
    for (int i = 0; i < reso; i++){
        for ( int j = 0; j < reso; j++){
            glVertex3f(x[j][i],y[j][i],z[j][i]);
            glVertex3f(x[j+1][i],y[j+1][i],z[j+1][i]);

            glVertex3f(x[j+1][i+1],y[j+1][i+1],z[j+1][i+1]);
            glVertex3f(x[j][i+1],y[j][i+1],z[j][i+1]);
        }
    }

    glEnd();
}

// define material color properties for front and back side
void SetMaterialColor( int side, float r, float g, float b){
    float	amb[4], dif[4], spe[4];
    int	i, mat;

    dif[0] = r; // diffuse color as defined by r,g, and b
    dif[1] = g;
    dif[2] = b;
    for( i=0; i<3; i++){
        amb[i] = .1 * dif[i]; // ambient color is 10 percent of diffuse
        spe[i] = .5; // specular color is just white / gray
    }
    amb[3] = dif[3] = spe[3] = 1.0; // alpha component is always 1
    switch( side){
    case 1:	mat = GL_FRONT; break;
    case 2:	mat = GL_BACK; break;
    default: mat = GL_FRONT_AND_BACK; break;
    }
    glMaterialfv( mat, GL_AMBIENT, amb); // define ambient, diffuse and specular components
    glMaterialfv( mat, GL_DIFFUSE, dif);
    glMaterialfv( mat, GL_SPECULAR, spe);
    glMaterialf( mat, GL_SHININESS, 50.0); // Phong constant for the size of highlights
}


OGLWidget::OGLWidget(QWidget *parent) // constructor
    : QOpenGLWidget(parent)
{
    // Setup the animation timer to fire every x msec
    animtimer = new QTimer(this);
    animtimer->start( 50 );

    // Everytime the timer fires, the animation is going one step forward
    connect(animtimer, SIGNAL(timeout()), this, SLOT(stepAnimation()));

    animstep = 0;
}

OGLWidget::~OGLWidget() // destructor
{
}

void OGLWidget::stepAnimation()
{
    animstep++;    // Increase animation steps
    update();      // Trigger redraw of scene with paintGL
}

void OGLWidget::initializeGL() // initializations to be called once
{
    initializeOpenGLFunctions();

    InitLightingAndProjection(); // define light sources and projection
}

void OGLWidget::paintGL() // draw everything, to be called repeatedly
{
    glEnable(GL_NORMALIZE); // this is necessary when using glScale (keep normals to unit length)

    // set background color
    glClearColor(0, .1, .1, .1); // bright blue
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // draw the scene
    glMatrixMode( GL_MODELVIEW);
    glLoadIdentity();				// Reset The Current Modelview Matrix




    //Torus

    glPushMatrix();
    glTranslated(.0,.0,.0);

    SetMaterialColor( 3, 0.0, 0.3, 1.0);
    glScaled( 1.2, 1.2, 1.2);
    glRotated(angleTorus1, 1, 1, -1);
    angleTorus1 += 1;
    DrawTorus(1,12);

    SetMaterialColor( 3, 0.1, 0.4, 0.8);
    glScaled( 1.0, 1., 1.);
    glRotated(angleTorus1, 1, 1, -1);
    //angleTorus2 += 1;
    DrawTorus(1,10);

    SetMaterialColor( 3, .2, .5, .6);
    glScaled( .8, .8, .8);
    glRotated(angleTorus1, 1, 1, -1);
    //angleTorus3 += 1;
    DrawTorus(1,10);

    SetMaterialColor( 3, .3, 0.6, 0.4);
    glScaled( .8, .8, .8);
    glRotated(angleTorus1, 1, 1, -1);
    //angleTorus3 += 1;
    DrawTorus(1,10);

    SetMaterialColor( 3, .4, .7, 0.2);
    glScaled( .8, .8, .8);
    glRotated(angleTorus1, 1, 1, -1);
    //angleTorus3 += 1;
    DrawTorus(1,10);



    //Mobius Strip

    SetMaterialColor( 3, 0.9, 0.2, 0.1);
    glScaled( 2.5,2.5,2.5);
    glRotated(beta, 1, 0, -1);
    beta+=3;
    DrawMoebius();


    //Sphere

    SetMaterialColor( 3, 1., .8, 0.1);
    glTranslated(3,0,0);
    glScaled( 0.4, 0.4, 0.4);
    //glRotated(alpha, 1, 1, -1);
    //alpha += 1;
    DrawSphere(1,0);


    //Cube

    glTranslated(-15,3,0);
    glScaled( .8, .8, .8);
    //glRotated(beta, 1, 1, 1);
    //beta += 3;
    DrawCube();

    //Pyramid

    glTranslated(10,5,0);
    glScaled( 1.0, 1.0, 1.0);
    //glRotated(beta, 0, 1, 1);
    DrawPyramid();


    // Cylinder

    glTranslated(0,-15,0);
    glScaled( 1.0, 1.0, 1.0);
    //glRotated(beta, 0, 1, 1);
    DrawCylinder();

    glPopMatrix();





    //Sphere
/*
    glPushMatrix();
    SetMaterialColor( 3, 1., .8, 0.1);
    glTranslated(4.5,0,0);
    glScaled( .8, .8, .8);
    DrawSphere(1,0);
    glRotated(alpha, alpha, alpha, -1 * alpha);
    alpha += 1;
    glPopMatrix();
*/



    //Cube
/*
    glPushMatrix();
    glTranslated(-2.0,-2.0,-2.0);
    glScaled( 1.5, 1.5, 1.5);
    DrawCube();
    glRotated(alpha, 0, 1, 1);
    glPopMatrix();
*/





    //Pyramid
/*
    glPushMatrix();
    glTranslated(2.0,2.0,2.0);
    glScaled( 2.0, 2.0, 2.0);
    DrawPyramid();
    glRotated(alpha, 0, 1, 1);
    glPopMatrix();
*/


    /*
    glTranslated( 0 ,0 ,-10.0);     // Move 10 units backwards in z, since camera is at origin
    glScaled( 1.0, 1.0, 1.0);       // scale objects
    glRotated( alpha, 0, 3, 1);     // continuous rotation
    alpha += 5; // Speed of rotation
*/
    // define color: 1=front, 2=back, 3=both, followed by r, g, and b
    SetMaterialColor( 1, 1.0, .2, .2);  // front color is red
    SetMaterialColor( 2, .2, .2, 1.0); // back color is blue

    //draw a cylinder with default resolution
    glPushMatrix();

    glPopMatrix();

    glRotated(-90,0,0,0);













    glFlush();
}

void OGLWidget::resizeGL(int w, int h) // called when window size is changed
{
    // adjust viewport transform
    glViewport(0,0,w,h);
}
