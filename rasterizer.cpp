#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "hw2_types.h"
#include "hw2_math_ops.h"
#include "hw2_file_ops.h"
#include <iostream>
#include <vector>

#define pi 3.14159265


Camera cameras[100];
int numberOfCameras = 0;

Model models[1000];
int numberOfModels = 0;

Color colors[100000];
int numberOfColors = 0;

Translation translations[1000];
int numberOfTranslations = 0;

Rotation rotations[1000];
int numberOfRotations = 0;

Scaling scalings[1000];
int numberOfScalings = 0;

Vec3 vertices[100000];
int numberOfVertices = 0;

Color backgroundColor;

// backface culling setting, default disabled
int backfaceCullingSetting = 0;

Color **image;

using namespace std;

/*
	Initializes image with background color
*/
void initializeImage(Camera cam) {
    int i, j;

    for (i = 0; i < cam.sizeX; i++)
        for (j = 0; j < cam.sizeY; j++) {
            image[i][j].r = backgroundColor.r;
            image[i][j].g = backgroundColor.g;
            image[i][j].b = backgroundColor.b;

        }
}

void color_assign(Color &eski,Color yeni) {
    eski.r = yeni.r;
    eski.g = yeni.g;
    eski.b = yeni.b;
}

void color_mult(Color &eski,double x) {
    eski.r *= x;
    eski.g *= x;
    eski.b *= x;
}

Color color_add(Color eski,Color yeni) {
    Color result;
    result.r = eski.r + yeni.r;
    result.g = eski.g + yeni.g;
    result.b = eski.b + yeni.b;
    return result;
}

Color color_sub(Color eski,Color yeni) {
    Color result;
    result.r = eski.r - yeni.r;
    result.g = eski.g - yeni.g;
    result.b = eski.b - yeni.b;
    return result;
}

double line_equation(Vec3 start, Vec3 end, int x, int y){
    return x*(start.y-end.y) + y*(end.x - start.x) + start.x*end.y - start.y*end.x;
}

double line_equation(Vec3 start, Vec3 end, double x, double y){
    return x*(start.y-end.y) + y*(end.x - start.x) + start.x*end.y - start.y*end.x;
}

int minindex(double x, double y, double z){
    if(abs(x) <= abs(y) && abs(x) < abs(z)) return 1;
    else if(abs(y) < abs(z) && abs(y) < abs(x)) return 2;
    else if(abs(z) <= abs(x) && abs(z) <= abs(y)) return 3;
}

/*
	Transformations, culling, rasterization are done here.
	You can define helper functions inside this file (rasterizer.cpp) only.
	Using types in "hw2_types.h" and functions in "hw2_math_ops.cpp" will speed you up while working.
*/
void forwardRenderingPipeline(Camera cam) {    
    for(int i = 0; i < numberOfModels; i++){
        for(int j = 0; j < models[i].numberOfTriangles; j++){
            Vec3 temp[3];
            for(int k = 0; k < 3; k++){
                Vec3 myvertice = vertices[ models[i].triangles[j].vertexIds[k] ];
                //cout << "Original: " << myvertice.x << " " << myvertice.y << " " << myvertice.z << endl;
                for(int l = 0; l < models[i].numberOfTransformations; l++){
                    if(models[i].transformationTypes[l] == 't'){
                        myvertice.x += translations[ models[i].transformationIDs[l] ].tx;
                        myvertice.y += translations[ models[i].transformationIDs[l] ].ty;
                        myvertice.z += translations[ models[i].transformationIDs[l] ].tz;
                    }

                    else if(models[i].transformationTypes[l] == 's'){
                        myvertice.x *= scalings[ models[i].transformationIDs[l] ].sx;
                        myvertice.y *= scalings[ models[i].transformationIDs[l] ].sy;
                        myvertice.z *= scalings[ models[i].transformationIDs[l] ].sz;                       
                    }

                    else{
                        double a = rotations[ models[i].transformationIDs[l] ].ux;
                        double b = rotations[ models[i].transformationIDs[l] ].uy;
                        double c = rotations[ models[i].transformationIDs[l] ].uz;
                        double angle = rotations[ models[i].transformationIDs[l] ].angle;
                        double ud = sqrt(b*b + c*c);
                        
                        double orgy = myvertice.y;
                        myvertice.y = c * myvertice.y/ud - myvertice.z * b/ud;
                        myvertice.z = c * myvertice.z/ud + orgy * b/ud;

                        double orgx = myvertice.x;
                        myvertice.x = ud * myvertice.x + a * myvertice.z;
                        myvertice.z = ud * myvertice.z - a * orgx;

                        orgx = myvertice.x;
                        myvertice.x = cos(angle * pi / 180) * myvertice.x - sin(angle * pi / 180) * myvertice.y;
                        myvertice.y = cos(angle * pi / 180) * myvertice.y + sin(angle * pi / 180) * orgx;

                        orgx = myvertice.x;
                        myvertice.x = ud * myvertice.x - a * myvertice.z;
                        myvertice.z = ud * myvertice.z + a * orgx;

                        orgy = myvertice.y;
                        myvertice.y = c * myvertice.y/ud + myvertice.z * b/ud;
                        myvertice.z = c * myvertice.z/ud - orgy * b/ud;
                    }
                    //cout << "Final: " << myvertice.x << " " << myvertice.y << " " << myvertice.z << endl;
                }
                double cam_matrix[4][4] = { {cam.u.x, cam.u.y, cam.u.z, -(cam.u.x * cam.pos.x + cam.u.y * cam.pos.y + cam.u.z * cam.pos.z)},
                                           {cam.v.x, cam.v.y, cam.v.z, -(cam.v.x * cam.pos.x + cam.v.y * cam.pos.y + cam.v.z * cam.pos.z)},
                                           {cam.w.x, cam.w.y, cam.w.z, -(cam.w.x * cam.pos.x + cam.w.y * cam.pos.y + cam.w.z * cam.pos.z)},
                                           {0.0,0.0,0.0,1.0} };
                double cam_perspective[4][4] =  { {2*cam.n/(cam.r - cam.l), 0, (cam.r + cam.l)/(cam.r - cam.l), 0},
                                                  {0, 2*cam.n/(cam.t - cam.b), (cam.t + cam.b)/(cam.t - cam.b), 0},
                                                  {0, 0, (cam.f + cam.n)/(cam.n - cam.f), (2*cam.f*cam.n)/(cam.n - cam.f)},
                                                  {0, 0, -1, 0} };
                double result_matrix[4][4] = {0};
                double vertice_matrix[4] = {myvertice.x, myvertice.y, myvertice.z, 1};
                double new_vertice[4] = {0};                                 
                multiplyMatrixWithMatrix(result_matrix, cam_perspective, cam_matrix);
                multiplyMatrixWithVec4d(new_vertice, result_matrix, vertice_matrix);
                

                temp[k].x = new_vertice[0]/new_vertice[3];  //perspective divide
                temp[k].y = new_vertice[1]/new_vertice[3];
                temp[k].z = new_vertice[2]/new_vertice[3];
                temp[k].colorId = myvertice.colorId;
            }
            Vec3 c_a = subtractVec3(temp[2], temp[0]);
            Vec3 b_a = subtractVec3(temp[1], temp[0]);
            Vec3 normal = normalizeVec3( crossProductVec3(c_a, b_a) );
            Vec3 viewing = normalizeVec3(temp[0]); 
            if(!backfaceCullingSetting || dotProductVec3(viewing, normal) <= 0){

                temp[0].x = ((temp[0].x + 1) * cam.sizeX - 1)/2;
                temp[0].y = ((temp[0].y + 1) * cam.sizeY - 1)/2;
                temp[1].x = ((temp[1].x + 1) * cam.sizeX - 1)/2;
                temp[1].y = ((temp[1].y + 1) * cam.sizeY - 1)/2;
                temp[2].x = ((temp[2].x + 1) * cam.sizeX - 1)/2;
                temp[2].y = ((temp[2].y + 1) * cam.sizeY - 1)/2;

                if(models[i].type == 1){
                    for(int hiter = 0; hiter < cam.sizeY; hiter++){
                        for(int witer = 0; witer < cam.sizeX; witer++){
                            double aa = line_equation(temp[2], temp[1], temp[0].x, temp[0].y);
                            double alpha = line_equation(temp[2], temp[1], witer, hiter) / aa;
                            double beta = line_equation(temp[0], temp[2], witer, hiter) / line_equation(temp[0], temp[2], temp[1].x, temp[1].y);
                            double gama = line_equation(temp[1], temp[0], witer, hiter) / line_equation(temp[1], temp[0], temp[2].x, temp[2].y);
                            if(alpha >= 0 && beta >= 0 && gama >= 0 && alpha <= 1 && beta <=1 && gama <= 1){
                                if(alpha > 10000000) cout << aa << endl;
                                image[witer][hiter].r = round(alpha * colors[temp[0].colorId].r + beta * colors[temp[1].colorId].r + gama * colors[temp[2].colorId].r);
                                image[witer][hiter].g = round(alpha * colors[temp[0].colorId].g + beta * colors[temp[1].colorId].g + gama * colors[temp[2].colorId].g);
                                image[witer][hiter].b = round(alpha * colors[temp[0].colorId].b + beta * colors[temp[1].colorId].b + gama * colors[temp[2].colorId].b);
                            }
                        }
                    }
                }
                else{
                    for(int liter = 0; liter < 3; liter++){
                        int x0 = round(temp[liter].x);
                        int x1 = round(temp[(liter+1)%3].x);
                        int y0 = round(temp[liter].y);
                        int y1 = round(temp[(liter+1)%3].y);
                        Color start_color;
                        Color end_color;
                        Color current_color;
                        Color dc;
                        double diff;
                        int tempp,hiter,witer;
                        if(x1 == x0) {  //infinite slope
                            if(y1 < y0) { //minus infinite
                                tempp = x0;
                                x0 = x1;
                                x1 = tempp;
                                tempp = y0;
                                y0 = y1;
                                y1 = y0;
                                color_assign(end_color,colors[temp[liter].colorId]);
                                color_assign(start_color,colors[temp[(liter+1)%3].colorId]);
                            }
                            else {
                                color_assign(start_color,colors[temp[liter].colorId]);
                                color_assign(end_color,colors[temp[(liter+1)%3].colorId]);
                            }
                            color_assign(current_color,start_color);
                            hiter = y0;
                            color_assign(dc,color_sub(end_color,start_color));
                            color_mult(dc,1/static_cast<double>(y1-y0));
                            while(hiter <= y1) {
                                color_assign(image[x0][hiter],current_color);
                                hiter++;
                                color_assign(current_color,color_add(current_color,dc));
                            }
                        }

                        else {
                            if(x1 < x0) { 
                                tempp = x0;
                                x0 = x1;
                                x1 = tempp;
                                tempp = y0;
                                y0 = y1;
                                y1 = tempp;
                                color_assign(end_color,colors[temp[liter].colorId]);
                                color_assign(start_color,colors[temp[(liter+1)%3].colorId]);
                            }
                            else {
                                color_assign(start_color,colors[temp[liter].colorId]);
                                color_assign(end_color,colors[temp[(liter+1)%3].colorId]);
                            }
                            color_assign(current_color,start_color);
                            double m = static_cast<double>(y1-y0)/(x1-x0);
                            if(m >= 0 && m <=1) { //Case 1
                                witer = x0;
                                hiter = y0;
                                color_assign(dc,color_sub(end_color,start_color));
                                color_mult(dc,1/static_cast<double>(x1-x0));
                                diff = 2*(y0-y1)+(x1-x0);
                                while(witer <= x1) {
                                    color_assign(image[witer][hiter],current_color);
                                    if(diff < 0) {
                                        hiter++;
                                        diff+= 2*(y0-y1) + 2*(x1-x0);
                                    }
                                    else {
                                        diff+= 2*(y0-y1);
                                    }
                                    color_assign(current_color,color_add(current_color,dc));
                                    witer++;
                                }
                            }

                            else if(m > 1) { //Case 2
                                witer = x0;
                                hiter = y0;
                                color_assign(dc,color_sub(end_color,start_color));
                                color_mult(dc,1/static_cast<double>(y1-y0));
                                diff = (y0-y1)+2*(x1-x0);
                                while(hiter <= y1) {
                                    color_assign(image[witer][hiter],current_color);
                                    if(diff < 0) {
                                        diff+= 2*(x1-x0);
                                    }
                                    else {
                                        witer++;
                                        diff+= 2*(y0-y1) + 2*(x1-x0);
                                    }
                                    color_assign(current_color,color_add(current_color,dc));
                                    hiter++;
                                }                                
                            }
                            else if(m >= -1) { //Case 3
                                witer = x0;
                                hiter = y0;
                                color_assign(dc,color_sub(end_color,start_color));
                                color_mult(dc,1/static_cast<double>(x1-x0));
                                diff = 2*(y0-y1)+(x0-x1);
                                while(witer <= x1) {
                                    color_assign(image[witer][hiter],current_color);
                                    if(diff < 0) {
                                        diff+= 2*(y0-y1);
                                    }
                                    else {
                                        hiter--;
                                        diff+= 2*(y0-y1) + 2*(x0-x1);
                                    }
                                    color_assign(current_color,color_add(current_color,dc));
                                    witer++;
                                }
                            }
                            else { //Case 4
                                witer = x0;
                                hiter = y0;
                                color_assign(dc,color_sub(end_color,start_color));
                                color_mult(dc,1/static_cast<double>(y0-y1));
                                diff = (y0-y1)+2*(x0-x1);
                                while(hiter >= y1) {
                                    color_assign(image[witer][hiter],current_color);
                                    if(diff < 0) {
                                        witer++;
                                        diff+= 2*(y0-y1) + 2*(x0-x1);
                                        
                                    }
                                    else {
                                        diff+= 2*(x0-x1);
                                    }
                                    color_assign(current_color,color_add(current_color,dc));
                                    hiter--;
                                }                                
                            }                            
                        }
                    }
                }
            }
        }
    }
}


int main(int argc, char **argv) {
    int i, j;

    if (argc < 2) {
        std::cout << "Usage: ./rasterizer <scene file> <camera file>" << std::endl;
        return 1;
    }

    // read camera and scene files
    readSceneFile(argv[1]);
    readCameraFile(argv[2]);

    image = 0;

    for (i = 0; i < numberOfCameras; i++) {

        // allocate memory for image
        if (image) {
			for (j = 0; j < cameras[i].sizeX; j++) {
		        delete image[j];
		    }

			delete[] image;
		}

        image = new Color*[cameras[i].sizeX];

        if (image == NULL) {
            std::cout << "ERROR: Cannot allocate memory for image." << std::endl;
            exit(1);
        }

        for (j = 0; j < cameras[i].sizeX; j++) {
            image[j] = new Color[cameras[i].sizeY];
            if (image[j] == NULL) {
                std::cout << "ERROR: Cannot allocate memory for image." << std::endl;
                exit(1);
            }
        }


        // initialize image with basic values
        initializeImage(cameras[i]);

        // do forward rendering pipeline operations
        forwardRenderingPipeline(cameras[i]);

        // generate PPM file
        writeImageToPPMFile(cameras[i]);

        // Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
        // Notice that os_type is not given as 1 (Ubuntu) or 2 (Windows), below call doesn't do conversion.
        // Change os_type to 1 or 2, after being sure that you have ImageMagick installed.
        convertPPMToPNG(cameras[i].outputFileName, 99);
    }

    return 0;

}
