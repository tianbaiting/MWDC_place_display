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
    double x, y, z,d;
    while (infile >> x >> y >> z>>d) {
        data.push_back({x, y, z,d});
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

    for (int i = 0; i < data.size(); ++i) {
    if (/* condition */data[i][3] == 3)
    {
        std::vector<double> coord = {data[i][0] - x, data[i][1] - y, data[i][2] - z};
        double dist;
        dist = distance(rotationMatrix, coord)[2];
        loss += (dist-d) * (dist-d);
        }
    else{
        std::vector<double> coord = {data[i][0] - x, data[i][1] - y, data[i][2] - z};
        double dist;
        dist = distance(rotationMatrix, coord)[data[i][3]];
        loss += dist * dist;
        }
    }



    return loss;
}

std::vector<double> chi2(double x,double y,double z,double d,double roll, double pitch, double yaw,std::vector<std::vector<double>> data){

    std::vector<double> chisqure;
    // 计算旋转矩阵
    std::vector<std::vector<double>> rotationMatrix = eulerToRotationMatrix(roll, pitch, yaw);


        for (int i = 0; i < data.size(); ++i) {
    if (/* condition */data[i][3] == 3)
    {
        std::vector<double> coord = {data[i][0] - x, data[i][1] - y, data[i][2] - z};
        double dist;
        dist = distance(rotationMatrix, coord)[2]-d;
        chisqure.push_back( dist * dist);
        }
    else{
        std::vector<double> coord = {data[i][0] - x, data[i][1] - y, data[i][2] - z};
        double dist;
        dist = distance(rotationMatrix, coord)[data[i][3]];
        chisqure.push_back( dist * dist);
        }
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

std::vector<double> BestPrameter(const std::string& filename) {
    // 读取 mwdc1.txt 文件中的数据
    globalData = readCoordinates(filename);

    // 创建 TMinuit 对象
    TMinuit minuit(7);
    minuit.SetFCN(fcn);

    //固定参数
    minuit.FixParameter(3);
  // 设置初始参数值和步长

//MWDC1
    double vstart[7] = {2000.0 , 
    602, 
    5000.0 ,
    50.0, 
    3.14/2 ,
    3.14- 25.0/180.0*3.14, 
    0.0};
    // MWDC0
    // double vstart[7] = {2000.0 , 
    // 768, 
    // 6000.0 ,
    // 100, 
    // 0.00411001 , 
    // 0.00753156, 
    // 0.0};

    // etof
    // double vstart[7] = {2000.0 , 766, 5000.0 ,187.0, 3.14 , 0.0, 1.5709};
    //MWDC



    double step[7] = {0.1, 0.1,0.1, 0.1, 0.1, 0.1, 0.1};


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

    return {x, y, z, d, roll, pitch, yaw};
}


// 反向旋转矩阵的列向量
void reverseRotationMatrixColumns(std::vector<std::vector<double>>& rotationMatrix, const std::vector<std::vector<double>>& globalData,double x,double y,double z) {
    for (int i = 0; i < 3; ++i) {
        double sum = 0.0;
        for (int j = 0; j < globalData.size(); ++j) {
            sum += 
              rotationMatrix[0][i] * (globalData[j][0]-x) 
            + rotationMatrix[1][i] * (globalData[j][1]-y) 
            + rotationMatrix[2][i] * (globalData[j][2]-z);
        }
        if (sum < 0) {
            for (int row = 0; row < 3; ++row) {
                rotationMatrix[0][i] = -rotationMatrix[0][i];
                rotationMatrix[1][i] = -rotationMatrix[1][i];
                rotationMatrix[2][i] = -rotationMatrix[2][i];
            }
            std::cout << "Column " << i << " was reversed." << std::endl;
        }
  
    }
}

int MWDC1_place() {

    // 读取 mwdc1.txt 文件中的数据

    std::string detectorName = "MWDC1";

    std::string filename1 = "./data_measure/" + detectorName + ".txt";
    std::vector<double> params = BestPrameter(filename1);

    double x = params[0];
    double y = params[1];
    double z = params[2];
    double d = params[3];
    double roll = params[4];
    double pitch = params[5];
    double yaw = params[6];


 std::vector<double> chi2_values = chi2(x,y,z,d,roll,pitch,yaw,globalData);

        // 创建 ROOT 文件
    std::string filename2 = "./root_chi2/" + detectorName + "_chi2_values.root";
    TFile *file = new TFile(filename2.c_str(), "RECREATE");

    // 创建 TTree
    TTree *tree = new TTree("chi2_tree", "Chi2 Values Tree");

    // 创建变量来存储 chi2 值
    double chi2_value;

    // 将变量添加到 TTree
    tree->Branch("chi2_value", &chi2_value, "chi2_value/D");

    // 填充 TTree
    for (double value : chi2_values) {
        chi2_value = value;
        std::cout << "chi2" <<detectorName<<" "  << chi2_value << std::endl;
        tree->Fill();
    }

    // 写入 ROOT 文件
    tree->Write();
    file->Close();

    
    // 生成旋转矩阵
    std::vector<std::vector<double>> rotationMatrix = eulerToRotationMatrix(roll, pitch, yaw);

// 反向旋转矩阵的列向量
    reverseRotationMatrixColumns(rotationMatrix, globalData,x,y,z);

    x= x + rotationMatrix[0][0] * 38.1 + rotationMatrix[0][1] * 138.1 + rotationMatrix[0][2] * 38.1;
    y= y + rotationMatrix[1][0] * 38.1 + rotationMatrix[1][1] * 138.1 + rotationMatrix[1][2] * 38.1;
    z= z + rotationMatrix[2][0] * 38.1 + rotationMatrix[2][1] * 138.1 + rotationMatrix[2][2] * 38.1;


    // 将 x, y, z 和旋转矩阵输出到文本文件
    std::ofstream outfile("output_" + detectorName+".txt");
    outfile << x << " " << y << " " << z << std::endl;
    for (const auto& row : rotationMatrix) {
        for (double value : row) {
            outfile << value << " ";
        }
        outfile << std::endl;
    }
    outfile.close();


   


    return 0;

}