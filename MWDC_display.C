#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TGLViewer.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TText.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <TColor.h>
#include <TMinuit.h>
#include <TMarker.h>
#include <TLatex.h>
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
    double x, y, z,d;
    while (infile >> x >> y >> z>>d) {
        data.push_back({x, y, z,d});
    }
    return data;
}

// 读取 chi2 值
std::vector<double> readChi2Values(const std::string& chi_filename) {
    std::vector<double> chi2_values;
    TFile *file = TFile::Open(chi_filename.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file " << chi_filename << std::endl;
        return chi2_values;
    }

    TTree *tree = nullptr;
    file->GetObject("chi2_tree", tree);
    if (!tree) {
        std::cerr << "Error getting TTree chi2_tree from file" << std::endl;
        file->Close();
        return chi2_values;
    }

    double chi2_value;
    tree->SetBranchAddress("chi2_value", &chi2_value);

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        chi2_values.push_back(chi2_value);
    }

    file->Close();
    return chi2_values;
}


void createSpheres(TGeoManager* geom, const std::vector<std::vector<double>>& coords, double radius, const std::string& chi_filename) {
    TGeoVolume *top = geom->GetTopVolume();
    TGeoMedium *medium = nullptr; // 使用默认介质

   // 读取 chi2 值
    std::vector<double> chi2_values = readChi2Values(chi_filename);

     // 定义颜色数组
    int colors[] = {kRed, kGreen, kYellow,kOrange};

    for (size_t i = 0; i < coords.size(); ++i) {
        TGeoVolume *sphere = geom->MakeSphere("SPHERE", medium, 0, radius);
        sphere->SetLineColor(coords[i][3] );
        TGeoTranslation *trans = new TGeoTranslation(coords[i][0], coords[i][1], coords[i][2]);
        top->AddNode(sphere, 0, trans);
                // 绘制误差
        if (i < chi2_values.size()) {
            double error = sqrt(chi2_values[i]);
            double textX = coords[i][0] + 100;
            double textY = coords[i][1] + 100;
            double textZ = coords[i][2];
            //  TText3D *text = new TText3D(textX, textY, textZ, Form("%.2f", error));
            // text->SetTextSize(0.02);
            // text->SetTextColor(kBlack);
            // text->Draw();

        }
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

          
   std::vector<std::string> detectorNames = {"eTOF", "MWDC0", "MWDC1", "T0", "target"};
//    std::vector<std::string> detectorNames = {"eTOF", "MWDC0", "MWDC1", "T0", "target"};



    for (const auto& detectorName : detectorNames) {
        std::string filename1 = "./data_measure/" + detectorName + ".txt";
        std::vector<std::vector<double>> coordinates = readCoordinates(filename1);
        std::string chi_filename = "./root_chi2/" + detectorName + "_chi2_values.root";
        createSpheres(geom, coordinates, 38.1, chi_filename);

        std::string directionFilename = "direction_" + detectorName + ".txt";
        std::vector<std::vector<double>> data = readDirectionTarget(directionFilename);
        createGeometry(geom, data, ("BOX_" + detectorName).c_str());
    }
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