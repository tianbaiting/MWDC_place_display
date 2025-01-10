#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TGLViewer.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <TColor.h>

// 读取文件并返回数据
std::vector<std::vector<double>> readDirectionTarget(const std::string& filename) {
    std::ifstream infile(filename);
    std::vector<std::vector<double>> data;
    double x, y, z;
    while (infile >> x >> y >> z) {
        data.push_back({x, y, z});
    }
    return data;
}

// 创建几何体并添加长方体
void createGeometry(TGeoManager* geom, const std::vector<std::vector<double>>& data, const char* boxName) {
    // 检查数据是否正确读取
    if (data.size() != 5) {
        std::cerr << "Error: direction_target.txt should contain exactly 5 lines of data." << std::endl;
        return;
    }

    // 获取顶点位置
    double vertex[3] = {data[0][0], data[0][1], data[0][2]};

    // 获取旋转矩阵
    double dir1[3] = {data[1][0], data[1][1],data[1][2]};
    double dir2[3] = {data[2][0], data[2][1],data[2][2]};   
    double dir3[3] = {data[3][0], data[3][1],data[3][2]};

    // 获取轴边长 半边长
    double lengths[3] = {data[4][0]/2, data[4][1]/2, data[4][2]/2};

    // 创建一个长方体
    TGeoVolume *box = geom->MakeBox(boxName, nullptr, lengths[0], lengths[1], lengths[2]);


   // 设置颜色和透明度
    box->SetLineColor(kBlue);
    box->SetTransparency(50); // 透明度范围为0-100，值越大越透明
    
    // double cx = vertex[0] + dir1[0] * lengths[0] + dir2[0] * lengths[1] + dir3[0] * lengths[2];
    // double cy = vertex[1] + dir1[1] * lengths[0] + dir2[1] * lengths[1] + dir3[1] * lengths[2];
    // double cz = vertex[2] + dir1[2] * lengths[0] + dir2[2] * lengths[1] + dir3[2] * lengths[2];
        
    double cx = vertex[0] + dir1[0] * lengths[0] + dir1[1] * lengths[1] + dir1[2] * lengths[2];
    double cy = vertex[1] + dir2[0] * lengths[0] + dir2[1] * lengths[1] + dir2[2] * lengths[2];
    double cz = vertex[2] + dir3[0] * lengths[0] + dir3[1] * lengths[1] + dir3[2] * lengths[2];
    // 创建变换矩阵
    TGeoTranslation *trans = new TGeoTranslation(cx, cy, cz);
    // TGeoTranslation *trans = new TGeoTranslation(0,0,0);
    TGeoRotation *rot = new TGeoRotation();
    Double_t rotationMatrix[9] = {dir1[0], dir1[1], dir1[2], dir2[0], dir2[1], dir2[2], dir3[0], dir3[1], dir3[2]};
    rot->SetMatrix(rotationMatrix);
    TGeoCombiTrans *comb = new TGeoCombiTrans(*trans, *rot);

    // 添加长方体到顶层体积
    geom->GetTopVolume()->AddNode(box, 0, comb);
}

// 读取文件并返回数据
std::vector<std::vector<double>> readCoordinates(const std::string& filename) {
    std::ifstream infile(filename);
    std::vector<std::vector<double>> data;
    double x, y, z;
    while (infile >> x >> y >> z) {
        data.push_back({x, y, z});
    }
    return data;
}

// 创建几何体并添加球体
void createSpheres(TGeoManager* geom, const std::vector<std::vector<double>>& coords, double radius) {
    TGeoVolume *top = geom->GetTopVolume();
    TGeoMedium *medium = nullptr; // 使用默认介质

     // 定义颜色数组
    int colors[] = {kRed, kGreen, kYellow};
    int numColors = sizeof(colors) / sizeof(colors[0]);

    for (size_t i = 0; i < coords.size(); ++i) {
        TGeoVolume *sphere = geom->MakeSphere("SPHERE", medium, 0, radius);
        sphere->SetLineColor(kRed); // 每四个点换一种颜色
        TGeoTranslation *trans = new TGeoTranslation(coords[i][0], coords[i][1], coords[i][2]);
        top->AddNode(sphere, 0, trans);
    }

}

// 创建坐标轴
void createAxes(TGeoManager* geom) {
    TGeoVolume *top = geom->GetTopVolume();
    TGeoMedium *medium = nullptr; // 使用默认介质

    // 创建X轴
    TGeoVolume *xAxis = geom->MakeTube("X_AXIS", medium, 0, 1, 4000);
    xAxis->SetLineColor(kRed);
    TGeoTranslation *xTrans = new TGeoTranslation(0, 0, 3000);
    top->AddNode(xAxis, 0,xTrans);

    // 创建Y轴
    TGeoVolume *yAxis = geom->MakeTube("Y_AXIS", medium, 0, 1, 1000);
    yAxis->SetLineColor(kGreen);
    TGeoRotation *yRot = new TGeoRotation();
    yRot->RotateX(90); // 将Y轴旋转90度，使其沿Y轴方向
    top->AddNode(yAxis, 0, yRot);

    // 创建Z轴
    TGeoVolume *zAxis = geom->MakeTube("Z_AXIS", medium, 0, 1, 1000);
    zAxis->SetLineColor(kBlue);
    TGeoRotation *zRot = new TGeoRotation();
    zRot->RotateY(90); // 将Z轴旋转90度，使其沿Z轴方向
    top->AddNode(zAxis, 0, zRot);
}

void MWDC_display(bool vis = true) {
    TGeoManager *geom = new TGeoManager("simple1", "Simple geometry");

    // 创建顶层体积
    TGeoVolume *top = geom->MakeBox("TOP", nullptr, 20000, 20000, 20000);
    geom->SetTopVolume(top);

          // 读取 data.txt 文件
    std::vector<std::vector<double>> coordinates = readCoordinates("./data_measure/mwdc0.txt");
    // 创建几何体并添加球体
    createSpheres(geom, coordinates, 25.2);
    // 读取 direction_MWDC1.txt 文件
    std::vector<std::vector<double>> data_mwdc0 = readDirectionTarget("direction_MWDC0.txt");
    // 创建几何体并添加第二个长方体
    createGeometry(geom, data_mwdc0, "BOX_MWDC0");

              
    std::vector<std::vector<double>> coordinates2 = readCoordinates("./data_measure/eTOF.txt");
    createSpheres(geom, coordinates2, 25.2);
    std::vector<std::vector<double>> data_eTOF= readDirectionTarget("direction_eTOF.txt");
    createGeometry(geom, data_eTOF, "BOX_eTOF");



    std::vector<std::vector<double>> coordinates3 = readCoordinates("./data_measure/MWDC1.txt");
    createSpheres(geom, coordinates3, 25.2);
    std::vector<std::vector<double>> data_MWDC1= readDirectionTarget("direction_MWDC1.txt");
    createGeometry(geom, data_MWDC1, "BOX_MWDC1");

    std::vector<std::vector<double>> coordinates4 = readCoordinates("./data_measure/target.txt");
    createSpheres(geom, coordinates4, 16.6);
    std::vector<std::vector<double>> data_target= readDirectionTarget("direction_target.txt");
    createGeometry(geom, data_target, "BOX_target");

    std::vector<std::vector<double>> coordinates5 = readCoordinates("./data_measure/T0.txt");
    createSpheres(geom, coordinates5, 16.6);
    std::vector<std::vector<double>> data_T0= readDirectionTarget("direction_T0.txt");
    createGeometry(geom, data_T0, "BOX_T0");
    // // 读取 direction_target.txt 文件
    // std::vector<std::vector<double>> data_target = readDirectionTarget("direction_target.txt");
    // // 创建几何体并添加第一个长方体
    // createGeometry(geom, data_target, "BOX_TARGET");

    

    // // 读取 direction_MWDC2.txt 文件
    // std::vector<std::vector<double>> data_mwdc2 = readDirectionTarget("direction_MWDC2.txt");
    // // 创建几何体并添加第三个长方体
    // createGeometry(geom, data_mwdc2, "BOX_MWDC2");

      // 创建坐标轴
    createAxes(geom);

    // 闭合几何体
    geom->CloseGeometry();

    // 绘制几何体
    geom->SetVisLevel(4);
    if (vis) {
        geom->GetTopVolume()->Draw("ogl");
    }

    // TGLViewer *viewer = dynamic_cast<TGLViewer*>(gPad->GetViewer3D());
    // viewer->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);
}

int main(int argc, char **argv) {
    TApplication app("ROOT Application", &argc, argv);
    MWDC_display();
    app.Run();
    return 0;
}