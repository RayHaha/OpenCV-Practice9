// Opencvtest.cpp : 定義主控台應用程式的進入點。
//

#include "stdafx.h"
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <iostream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <math.h>

using namespace cv;
using namespace std;

void Robert(Mat originimage, Mat resultimage, double threshhold);
void Prewitt(Mat originimage, Mat resultimage, double threshold);
void Sobel(Mat originimage, Mat resultimage, double threshold);
void FreiAndChen(Mat originimage, Mat resultimage, double threshold);
void Kirsch(Mat originimage, Mat resultimage, double threshold);
void Robinson(Mat originimage, Mat resultimage, double threshold);
void NevatiaAndBabu(Mat originimage, Mat resultimage, double threshold);

int _tmain(int argc, _TCHAR* argv[])
{
	Mat image = imread("lena.bmp");
	Mat grayimage;
	int pixel = 0;
	cvtColor(image,grayimage, CV_BGR2GRAY);
	
	Mat robertimage = grayimage.clone();
	Robert(grayimage, robertimage, 15);
	imwrite("robert.bmp", robertimage);

	Mat prewittimage = grayimage.clone();
	Prewitt(grayimage, prewittimage, 45);
	imwrite("prewitt.bmp", prewittimage);

	Mat sobelimage = grayimage.clone();
	Sobel(grayimage, sobelimage, 55);
	imwrite("sobel.bmp", sobelimage);

	Mat freiandchenimage = grayimage.clone();
	FreiAndChen(grayimage, freiandchenimage, 50);
	imwrite("freiandchen.bmp", freiandchenimage);

	Mat kirschimage = grayimage.clone();
	Kirsch(grayimage, kirschimage, 145);
	imwrite("kirsch.bmp", kirschimage);

	Mat robinsonimage = grayimage.clone();
	Robinson(grayimage, robinsonimage, 55);
	imwrite("robinson.bmp", robinsonimage);

	
	Mat nevatiaandbabuimage = grayimage.clone();
	NevatiaAndBabu(grayimage, nevatiaandbabuimage, 12900);
	imwrite("nevatiaandbabu.bmp", nevatiaandbabuimage);
	

	waitKey(0);
	return 0;
}

void Robert(Mat originimage, Mat resultimage, double threshold){
	for(int i=0; i<originimage.rows-1; i++){
		for(int j=0; j<originimage.cols-1; j++){
			int element00 = originimage.at<uchar>(i,j);
			int element01 = originimage.at<uchar>(i,j+1);
			int element10 = originimage.at<uchar>(i+1,j);
			int element11 = originimage.at<uchar>(i+1,j+1);
			int pixel1 = element00*(-1)+element11;
			int pixel2 = element01*(-1)+element10;
			double pixel = sqrt((pixel1*pixel1)+(pixel2*pixel2));
			if(pixel >= threshold){
				resultimage.at<uchar>(i,j) = 0;
			}else{
				resultimage.at<uchar>(i,j) = 255;
			}
		}
	}
}

void Prewitt(Mat originimage, Mat resultimage, double threshold){
	for(int i=1; i<originimage.rows-1; i++){
		for(int j=1; j<originimage.cols-1; j++){
			int p1 = 0;
			for(int k=j-1; k<j+2; k++){
				p1 = p1 - originimage.at<uchar>(i-1,k);
			}
			for(int k=j-1; k<j+2; k++){
				p1 = p1 + originimage.at<uchar>(i+1,k);
			}
			int p2 = 0;
			for(int k=i-1; k<i+2; k++){
				p2 = p2 - originimage.at<uchar>(k,j-1);
			}
			for(int k=i-1; k<i+2; k++){
				p2 = p2 + originimage.at<uchar>(k,j+1);
			}
			double pixel = sqrt((p1*p1)+(p2*p2));
			if(pixel >= threshold){
				resultimage.at<uchar>(i,j) = 0;
			}else{
				resultimage.at<uchar>(i,j) = 255;
			}
		}
	}
}

void Sobel(Mat originimage, Mat resultimage, double threshold){
	for(int i=1; i<originimage.rows-1; i++){
		for(int j=1; j<originimage.cols-1; j++){
			int p1 = 0;
			for(int k=j-1; k<j+2; k++){
				p1 = p1 - originimage.at<uchar>(i-1,k);
			}
			p1 = p1 - originimage.at<uchar>(i-1,j);
			for(int k=j-1; k<j+2; k++){
				p1 = p1 + originimage.at<uchar>(i+1,k);
			}
			p1 = p1 + originimage.at<uchar>(i+1,j);
			int p2 = 0;
			for(int k=i-1; k<i+2; k++){
				p2 = p2 - originimage.at<uchar>(k,j-1);
			}
			p2 = p2 - originimage.at<uchar>(i,j-1);
			for(int k=i-1; k<i+2; k++){
				p2 = p2 + originimage.at<uchar>(k,j+1);
			}
			p2 = p2 + originimage.at<uchar>(i,j+1);
			double pixel = sqrt((p1*p1)+(p2*p2));
			if(pixel >= threshold){
				resultimage.at<uchar>(i,j) = 0;
			}else{
				resultimage.at<uchar>(i,j) = 255;
			}
		}
	}
}

void FreiAndChen(Mat originimage, Mat resultimage, double threshold){
	for(int i=1; i<originimage.rows-1; i++){
		for(int j=1; j<originimage.cols-1; j++){
			double p1 = 0;
			p1 = (-1)*originimage.at<uchar>(i-1,j-1) + (-1)*sqrt(2)*originimage.at<uchar>(i-1,j) + (-1)*originimage.at<uchar>(i-1,j+1);
			p1 = p1 + originimage.at<uchar>(i+1,j-1) + sqrt(2)*originimage.at<uchar>(i+1,j) + originimage.at<uchar>(i+1,j+1);
			double p2 = 0;
			p2 = (-1)*originimage.at<uchar>(i-1,j-1) + (-1)*sqrt(2)*originimage.at<uchar>(i,j-1) + (-1)*originimage.at<uchar>(i+1,j-1);
			p2 = p2 + originimage.at<uchar>(i-1,j+1) + sqrt(2)*originimage.at<uchar>(i,j+1) + originimage.at<uchar>(i+1,j+1);
			double pixel = sqrt((p1*p1)+(p2*p2));
			if(pixel >= threshold){
				resultimage.at<uchar>(i,j) = 0;
			}else{
				resultimage.at<uchar>(i,j) = 255;
			}
		}
	}
}

void Kirsch(Mat originimage, Mat resultimage, double threshold){
	for(int i=1; i<originimage.rows-1; i++){
		for(int j=1; j<originimage.cols-1; j++){
			double k[8] = {0,0,0,0,0,0,0,0};
			for(int u=i-1; u<i+2; u++){
				for(int v=j-1; v<j+2; v++){
					if(u==i && v==j){
						
					}else{
						if(v==j+1){
							k[0] = k[0] + 5*originimage.at<uchar>(u,v);
						}else{
							k[0] = k[0] - 3*originimage.at<uchar>(u,v);
						}
					}
					
				}
			}
			for(int u=i-1; u<i+2; u++){
				for(int v=j-1; v<j+2; v++){
					if(u==i && v==j){
						
					}else{
						if((v==j+1 && u!=i+1) || (u==i-1 && v==j)){
							k[1] = k[1] + 5*originimage.at<uchar>(u,v);
						}else{
							k[1] = k[1] - 3*originimage.at<uchar>(u,v);
						}
					}
					
				}
			}
			for(int u=i-1; u<i+2; u++){
				for(int v=j-1; v<j+2; v++){
					if(u==i && v==j){
						
					}else{
						if(u==i-1){
							k[2] = k[2] + 5*originimage.at<uchar>(u,v);
						}else{
							k[2] = k[2] - 3*originimage.at<uchar>(u,v);
						}
					}
				}
			}
			for(int u=i-1; u<i+2; u++){
				for(int v=j-1; v<j+2; v++){
					if(u==i && v==j){
						
					}else{
						if((u==i-1 && v!=j+1) || (u==i && v==j-1)){
							k[3] = k[3] + 5*originimage.at<uchar>(u,v);
						}else{
							k[3] = k[3] - 3*originimage.at<uchar>(u,v);
						}
					}
				}
			}
			for(int u=i-1; u<i+2; u++){
				for(int v=j-1; v<j+2; v++){
					if(u==i && v==j){
						
					}else{
						if(v==j-1){
							k[4] = k[4] + 5*originimage.at<uchar>(u,v);
						}else{
							k[4] = k[4] - 3*originimage.at<uchar>(u,v);
						}
					}
				}
			}
			for(int u=i-1; u<i+2; u++){
				for(int v=j-1; v<j+2; v++){
					if(u==i && v==j){
						
					}else{
						if((v==j-1 && u!=i-1) || (v==j && u==i+1)){
							k[5] = k[5] + 5*originimage.at<uchar>(u,v);
						}else{
							k[5] = k[5] - 3*originimage.at<uchar>(u,v);
						}
					}
				}
			}
			for(int u=i-1; u<i+2; u++){
				for(int v=j-1; v<j+2; v++){
					if(u==i && v==j){
						
					}else{
						if(u==i+1){
							k[6] = k[6] + 5*originimage.at<uchar>(u,v);
						}else{
							k[6] = k[6] - 3*originimage.at<uchar>(u,v);
						}
					}
				}
			}
			for(int u=i-1; u<i+2; u++){
				for(int v=j-1; v<j+2; v++){
					if(u==i && v==j){
						
					}else{
						if((u==i+1 && v!=j-1) || (u==i && v==j+1)){
							k[7] = k[7] + 5*originimage.at<uchar>(u,v);
						}else{
							k[7] = k[7] - 3*originimage.at<uchar>(u,v);
						}
					}
				}
			}

			double largest = -100;
			for(int u=0; u<8; u++){
				if(k[u]>largest){
					largest = k[u];
				}
			}
			if(largest >= threshold){
				resultimage.at<uchar>(i,j) = 0;
			}else{
				resultimage.at<uchar>(i,j) = 255;
			}
		}
	}
}

void Robinson(Mat originimage, Mat resultimage, double threshold){
	for(int i=1; i<originimage.rows-1; i++){
		for(int j=1; j<originimage.cols-1; j++){
			double r[8] = {0,0,0,0,0,0,0,0};
			
			r[0] = (-1)*originimage.at<uchar>(i-1,j-1)+(-2)*originimage.at<uchar>(i,j-1)+(-1)*originimage.at<uchar>(i+1,j-1);
			r[0] = r[0] + originimage.at<uchar>(i-1,j+1)+2*originimage.at<uchar>(i,j+1)+originimage.at<uchar>(i+1,j+1);

			r[1] = (-1)*originimage.at<uchar>(i,j-1)+(-2)*originimage.at<uchar>(i+1,j-1)+(-1)*originimage.at<uchar>(i+1,j);
			r[1] = r[1] + originimage.at<uchar>(i,j+1)+2*originimage.at<uchar>(i+1,j+1)+originimage.at<uchar>(i-1,j);

			r[2] = (-1)*originimage.at<uchar>(i+1,j-1)+(-2)*originimage.at<uchar>(i+1,j)+(-1)*originimage.at<uchar>(i+1,j+1);
			r[2] = r[2] + originimage.at<uchar>(i-1,j-1)+2*originimage.at<uchar>(i-1,j)+originimage.at<uchar>(i-1,j+1);
			
			r[3] = (-1)*originimage.at<uchar>(i+1,j)+(-2)*originimage.at<uchar>(i+1,j+1)+(-1)*originimage.at<uchar>(i,j+1);
			r[3] = r[3] + originimage.at<uchar>(i,j-1)+2*originimage.at<uchar>(i-1,j-1)+originimage.at<uchar>(i-1,j);

			r[4] = 1*originimage.at<uchar>(i-1,j-1)+2*originimage.at<uchar>(i,j-1)+1*originimage.at<uchar>(i+1,j-1);
			r[4] = r[4] + (-1)*originimage.at<uchar>(i-1,j+1)+(-2)*originimage.at<uchar>(i,j+1)+(-1)*originimage.at<uchar>(i+1,j+1);

			r[5] = originimage.at<uchar>(i,j-1)+2*originimage.at<uchar>(i+1,j-1)+originimage.at<uchar>(i+1,j);
			r[5] = r[5] + (-1)*originimage.at<uchar>(i,j+1)+(-2)*originimage.at<uchar>(i+1,j+1)+(-1)*originimage.at<uchar>(i-1,j);

			r[6] = originimage.at<uchar>(i+1,j-1)+2*originimage.at<uchar>(i+1,j)+originimage.at<uchar>(i+1,j+1);
			r[6] = r[6] + (-1)*originimage.at<uchar>(i-1,j-1)+(-2)*originimage.at<uchar>(i-1,j)+(-1)*originimage.at<uchar>(i-1,j+1);
			
			r[7] = originimage.at<uchar>(i+1,j)+2*originimage.at<uchar>(i+1,j+1)+originimage.at<uchar>(i,j+1);
			r[7] = r[7] + (-1)*originimage.at<uchar>(i,j-1)+(-2)*originimage.at<uchar>(i-1,j-1)+(-1)*originimage.at<uchar>(i-1,j);


			double largest = -100;
			for(int u=0; u<8; u++){
				if(r[u]>largest){
					largest = r[u];
				}
			}
			if(largest >= threshold){
				resultimage.at<uchar>(i,j) = 0;
			}else{
				resultimage.at<uchar>(i,j) = 255;
			}
		}
	}
}

void NevatiaAndBabu(Mat originimage, Mat resultimage, double threshold){
	for(int i=2; i<originimage.rows-2; i++){
		for(int j=2; j<originimage.cols-2; j++){
			double r[6] = {0,0,0,0,0,0};

			// r[0] degree 0
			for(int u=i-2; u<i; u++){
				for(int v=j-2; v<j+3; v++){
					r[0] = r[0] + originimage.at<uchar>(u,v);
				}
			}
			for(int u=i+1; u<i+3; u++){
				for(int v=j-2; v<j+3; v++){
					r[0] = r[0] - originimage.at<uchar>(u,v);
				}
			}

			// r[1] degree 30
			for(int v=j-2; v<j+3; v++){
				r[1] = r[1] + 100*originimage.at<uchar>(i-2,v);
			}
			r[1] = r[1] + 100*originimage.at<uchar>(i-1,j-2)+100*originimage.at<uchar>(i-1,j-1)+100*originimage.at<uchar>(i-1,j)
				+78*originimage.at<uchar>(i-1,j+1)-32*originimage.at<uchar>(i-1,j+2);
			r[1] = r[1] + 100*originimage.at<uchar>(i,j-2)+92*originimage.at<uchar>(i,j-1)
				-92*originimage.at<uchar>(i,j+1)-100*originimage.at<uchar>(i,j+2);
			r[1] = r[1] + 32*originimage.at<uchar>(i+1,j-2)-78*originimage.at<uchar>(i+1,j-1)-100*originimage.at<uchar>(i+1,j)
				-100*originimage.at<uchar>(i+1,j+1)-100*originimage.at<uchar>(i+1,j+2);
			for(int v=j-2; v<j+3; v++){
				r[1] = r[1] - 100*originimage.at<uchar>(i+2,v);
			}

			// r[2] degree 60
			for(int v=j-2; v<j+1; v++){
				for(int u=i-2; u<i-2*(v-j)-1; u++){
					r[2] = r[2] + 100*originimage.at<uchar>(u,v);
				}
			}
			r[2] = r[2] + 32*originimage.at<uchar>(i-2,j+1)+92*originimage.at<uchar>(i-1,j)-78*originimage.at<uchar>(i-1,j+1)
				+78*originimage.at<uchar>(i+1,j-1)-92*originimage.at<uchar>(i+1,j)-32*originimage.at<uchar>(i+2,j-1);
			for(int v=j+2; v>j-1; v--){
				for(int u=i+2; u>i-2*(v-j)+1; u--){
					r[2] = r[2] - 100*originimage.at<uchar>(u,v);
				}
			}

			// r[3] degree -90
			for(int v=j-2; v<j; v++){
				for(int u=i-2; u<i+3; u++){
					r[3] = r[3] - 100*originimage.at<uchar>(u,v);
				}
			}
			for(int v=j+1; v<j+3; v++){
				for(int u=i-2; u<i+3; u++){
					r[3] = r[3] + 100*originimage.at<uchar>(u,v);
				}
			}

			// r[4] degree -60
			for(int v=j-2; v<j+1; v++){
				for(int u=i+2; u>i+2*(v-j)+1; u--){
					r[4] = r[4] - 100*originimage.at<uchar>(u,v);
				}
			}
			r[4] = r[4] + 32*originimage.at<uchar>(i-2,j-1)+92*originimage.at<uchar>(i-1,j)-78*originimage.at<uchar>(i-1,j-1)
				+78*originimage.at<uchar>(i+1,j+1)-92*originimage.at<uchar>(i+1,j)-32*originimage.at<uchar>(i+2,j+1);
			for(int v=j+2; v>j-1; v--){
				for(int u=i-2; u<i+2*(v-j)-1; u++){
					r[4] = r[4] + 100*originimage.at<uchar>(u,v);
				}
			}

			// r[5] degree -30
			for(int v=j-2; v<j+3; v++){
				r[1] = r[1] + 100*originimage.at<uchar>(i-2,v);
			}
			r[1] = r[1] + 100*originimage.at<uchar>(i-1,j+2)+100*originimage.at<uchar>(i-1,j+1)+100*originimage.at<uchar>(i-1,j)
				+78*originimage.at<uchar>(i-1,j-1)-32*originimage.at<uchar>(i-1,j-2);
			r[1] = r[1] + 100*originimage.at<uchar>(i,j+2)+92*originimage.at<uchar>(i,j+1)
				-92*originimage.at<uchar>(i,j-1)-100*originimage.at<uchar>(i,j-2);
			r[1] = r[1] + 32*originimage.at<uchar>(i+1,j+2)-78*originimage.at<uchar>(i+1,j+1)-100*originimage.at<uchar>(i+1,j)
				-100*originimage.at<uchar>(i+1,j-1)-100*originimage.at<uchar>(i+1,j-2);
			for(int v=j-2; v<j+3; v++){
				r[1] = r[1] - 100*originimage.at<uchar>(i+2,v);
			}



			double largest = -100;
			for(int u=0; u<6; u++){
				if(r[u]>largest){
					largest = r[u];
				}
			}
			if(largest >= threshold){
				resultimage.at<uchar>(i,j) = 0;
			}else{
				resultimage.at<uchar>(i,j) = 255;
			}
		}
	}
}