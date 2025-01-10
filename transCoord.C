#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TGLViewer.h>
#include <TH1F.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <TColor.h>
#include <TMinuit.h>
#include <TFile.h>
#include <TTree.h>

std::vector<std::vector<double>> readDirectionTarget(const std::string& filename) {
    std::ifstream infile(filename);
    std::vector<std::vector<double>> data1;
    double x, y, z;
    while (infile >> x >> y >> z) {
        data1.push_back({x, y, z});
    }
    return data1;
}

std::vector<double> trans_MWDC_to_Local(std::string detectorName,std::vector<double> coord)
{   
    std::string directionFilename = "direction_" + detectorName + ".txt";
    const std::vector<std::vector<double>> data = readDirectionTarget(directionFilename);
   // 获取顶点位置
    double vertex[3] = {data[0][0], data[0][1], data[0][2]};

    // 获取旋转矩阵
    double dir1[3] = {data[1][0], data[1][1],data[1][2]};
    double dir2[3] = {data[2][0], data[2][1],data[2][2]};   
    double dir3[3] = {data[3][0], data[3][1],data[3][2]};

    double cx = vertex[0] + dir1[0] * coord[0] + dir1[1] * coord[1] + dir1[2] * coord[2];
    double cy = vertex[1] + dir2[0] * coord[0] + dir2[1] * coord[1] + dir2[2] * coord[2];
    double cz = vertex[2] + dir3[0] * coord[0] + dir3[1] * coord[1] + dir3[2] * coord[2];
    return {cx,cy,cz};
}


std::vector<double> transCoord()
{
    std::string detectorName ="MWDC1";
    std::vector<double> coord = {0,0,0};
    std::vector<double> coord_local = trans_MWDC_to_Local(detectorName,coord);
    
    return coord_local;
}