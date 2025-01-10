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


std::vector<std::vector<double>> globalData;
// 读取文件并返回数据
std::vector<std::vector<double>>readCoordinates(const std::string& filename) {
    std::ifstream infile(filename);
    std::vector<std::vector<double>> data;
    double x, y, z;
    while (infile >> x >> y >> z) {
        data.push_back({x, y, z});
    }
    return data;
}

// 将欧拉角转换为旋转矩阵
std::vector<std::vector<double>> eulerToRotationMatrix(double roll, double pitch, double yaw) {
    // 计算各个角度的正弦和余弦值
    double c1 = cos(yaw);
    double s1 = sin(yaw);
    double c2 = cos(pitch);
    double s2 = sin(pitch);
    double c3 = cos(roll);
    double s3 = sin(roll);

    // 生成旋转矩阵
    std::vector<std::vector<double>> rotationMatrix(3, std::vector<double>(3));
    rotationMatrix[0][0] = c1 * c2;
    rotationMatrix[0][1] = c1 * s2 * s3 - s1 * c3;
    rotationMatrix[0][2] = c1 * s2 * c3 + s1 * s3;
    rotationMatrix[1][0] = s1 * c2;
    rotationMatrix[1][1] = s1 * s2 * s3 + c1 * c3;
    rotationMatrix[1][2] = s1 * s2 * c3 - c1 * s3;
    rotationMatrix[2][0] = -s2;
    rotationMatrix[2][1] = c2 * s3;
    rotationMatrix[2][2] = c2 * c3;

    return rotationMatrix;
}



 // 计算旋转后的单位矢量与给定坐标的内积
std::vector<double> distance(const std::vector<std::vector<double>>& rotationMatrix, const std::vector<double>& coord) {
    if (rotationMatrix.size() != 3 || rotationMatrix[0].size() != 3 || coord.size() != 3) {
        std::cerr << "Invalid input dimensions." << std::endl;
        return {0,0,0};
    }

    // 计算旋转后的坐标
    double x = rotationMatrix[0][0] * coord[0] + rotationMatrix[1][0] * coord[1] + rotationMatrix[2][0] * coord[2];
    double y = rotationMatrix[0][1] * coord[0] + rotationMatrix[1][1] * coord[1] + rotationMatrix[2][1] * coord[2];
    double z = rotationMatrix[0][2] * coord[0] + rotationMatrix[1][2] * coord[1] + rotationMatrix[2][2] * coord[2];

    // 返回旋转后的坐标数组
    return {x, y, z};
}

double LossFunction(double x,double y,double z,double d,double roll, double pitch, double yaw,std::vector<std::vector<double>> data){
    // 计算旋转矩阵
    std::vector<std::vector<double>> rotationMatrix = eulerToRotationMatrix(roll, pitch, yaw);

    double loss = 0.0;

    // 前四个坐标计算 distance(rotationMatrix, data[i] - {x, y, z})[0]^2
    for (int i = 0; i < 4; ++i) {
        std::vector<double> coord = {data[i][0] - x, data[i][1] - y, data[i][2] - z};
        double dist = distance(rotationMatrix, coord)[0];
        loss += dist * dist;
    }

    // 中间四个坐标计算 distance(rotationMatrix, data[i] - {x, y, z})[1]^2
    for (int i = 4; i < 8; ++i) {
        std::vector<double> coord = {data[i][0] - x, data[i][1] - y, data[i][2] - z};
        double dist = distance(rotationMatrix, coord)[1];
        loss += dist * dist;
    }

    // 后面四个坐标计算 distance(rotationMatrix, data[i] - {x, y, z})[2]^2
    for (int i = 8; i < 12; ++i) {
        std::vector<double> coord = {data[i][0] - x, data[i][1] - y, data[i][2] - z};
        double dist = distance(rotationMatrix, coord)[2];
        loss += dist * dist;
    }

    for(int i = 12 ; i<16; ++i){
    
        std::vector<double> coord = {data[i][0] - x, data[i][1] - y, data[i][2] - z};
        double dist = distance(rotationMatrix, coord)[2]-d;
        loss += dist * dist;
    }

    return loss;
}

std::vector<double> chi2(double x,double y,double z,double d,double roll, double pitch, double yaw,std::vector<std::vector<double>> data){

    std::vector<double> chisqure;
    // 计算旋转矩阵
    std::vector<std::vector<double>> rotationMatrix = eulerToRotationMatrix(roll, pitch, yaw);

    // 前四个坐标计算 distance(rotationMatrix, data[i] - {x, y, z})[0]^2
    for (int i = 0; i < 4; ++i) {
        std::vector<double> coord = {data[i][0] - x, data[i][1] - y, data[i][2] - z};
        double dist = distance(rotationMatrix, coord)[0];
        chisqure.push_back( dist * dist);

    }

    // 中间四个坐标计算 distance(rotationMatrix, data[i] - {x, y, z})[1]^2
    for (int i = 4; i < 8; ++i) {
        std::vector<double> coord = {data[i][0] - x, data[i][1] - y, data[i][2] - z};
        double dist = distance(rotationMatrix, coord)[1];
        chisqure.push_back( dist * dist);
        
    }

    // 后面四个坐标计算 distance(rotationMatrix, data[i] - {x, y, z})[2]^2
    for (int i = 8; i < 12; ++i) {
        std::vector<double> coord = {data[i][0] - x, data[i][1] - y, data[i][2] - z};
        double dist = distance(rotationMatrix, coord)[2];
        chisqure.push_back( dist * dist);

    }

    for(int i = 12 ; i<16; ++i){
    
        std::vector<double> coord = {data[i][0] - x, data[i][1] - y, data[i][2] - z};
        double dist = distance(rotationMatrix, coord)[2]-d;
        chisqure.push_back( dist * dist);
    }
    return chisqure;

}

// TMinuit 计算 LossFunction 的包装函数
void fcn(int& npar, double* deriv, double& f, double par[], int flag) {
    double x = par[0];
    double y = par[1];
    double z = par[2];
    double d = par[3];
    double roll = par[4];
    double pitch = par[5];
    double yaw = par[6];

    f = LossFunction(x, y, z,d, roll, pitch, yaw, globalData);
}



int MWDC1_place() {

  // 读取 mwdc1.txt 文件中的数据
    globalData = readCoordinates("./data_measure/MWDC1.txt");

    // 创建 TMinuit 对象
    TMinuit minuit(7);
    minuit.SetFCN(fcn);

    // 设置初始参数值和步长
    double vstart[7] = {2000.0 , 602, 5000.0 ,50.0, 3.14/2 ,3.14- 25.0/180.0*3.14, 0.0};
    double step[7] = {0.1, 0.1,0.1, 0.1, 0.1, 0.1,0.1};

    minuit.DefineParameter(0, "x", vstart[0], step[0], 0 , 0);
    minuit.DefineParameter(1, "y", vstart[1], step[1], 0 ,0);
    minuit.DefineParameter(2, "z", vstart[2], step[2],0,0);
    minuit.DefineParameter(3, "d", vstart[3], step[3],0,0);
    minuit.DefineParameter(4, "roll", vstart[4], step[4], 0, 2*M_PI);
    minuit.DefineParameter(5, "pitch", vstart[5], step[5], -M_PI/2, M_PI/2);
    minuit.DefineParameter(6, "yaw", vstart[6], step[6], 0,2* M_PI);

    // 执行最小化
    minuit.Migrad();

    // 获取最优参数值
    double x, y, z,d, roll, pitch, yaw;
    minuit.GetParameter(0, x, step[0]);
    minuit.GetParameter(1, y, step[1]);
    minuit.GetParameter(2, z, step[2]);
    minuit.GetParameter(3, d, step[3]);
    minuit.GetParameter(4, roll, step[4]);
    minuit.GetParameter(5, pitch, step[5]);
    minuit.GetParameter(6, yaw, step[6]);

    // 输出最优参数值
    std::cout << "Optimal parameters:" << std::endl;
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    std::cout << "z = " << z << std::endl;
    std::cout << "d = " << d << std::endl;
    std::cout << "roll = " << roll << std::endl;
    std::cout << "pitch = " << pitch << std::endl;
    std::cout << "yaw = " << yaw << std::endl;

    
    // 生成旋转矩阵
    std::vector<std::vector<double>> rotationMatrix = eulerToRotationMatrix(roll, pitch, yaw);

    // 将 x, y, z 和旋转矩阵输出到文本文件
    std::ofstream outfile("output.txt");
    outfile << x << " " << y << " " << z << std::endl;
    for (const auto& row : rotationMatrix) {
        for (double value : row) {
            outfile << value << " ";
        }
        outfile << std::endl;
    }
    outfile.close();


    std::vector<double> chi2_values = chi2(x,y,z,d,roll,pitch,yaw,globalData);

        // 创建 ROOT 文件
    TFile *file = new TFile("./root_chi2/MWDC1_chi2_values.root", "RECREATE");

    // 创建 TTree
    TTree *tree = new TTree("chi2_tree", "Chi2 Values Tree");

    // 创建变量来存储 chi2 值
    double chi2_value;

    // 将变量添加到 TTree
    tree->Branch("chi2_value", &chi2_value, "chi2_value/D");

    // 填充 TTree
    for (double value : chi2_values) {
        chi2_value = value;
        tree->Fill();
    }

    // 写入 ROOT 文件
    tree->Write();
    file->Close();


    return 0;

}