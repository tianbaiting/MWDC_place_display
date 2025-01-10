#include <TFile.h>
#include <TTree.h>
#include <iostream>

int tree_out() {
    // 打开 ROOT 文件
    TFile *file = TFile::Open("chi2_values.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file chi2_values.root" << std::endl;
        return -1;
    }

    // 获取 TTree
    TTree *tree = nullptr;
    file->GetObject("chi2_tree", tree);
    if (!tree) {
        std::cerr << "Error getting TTree chi2_tree from file" << std::endl;
        file->Close();
        return -1;
    }

    // 设置分支地址
    double chi2_value;
    tree->SetBranchAddress("chi2_value", &chi2_value);

    // 输出 TTree 数据
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        std::cout << "Entry " << i << ": chi2_value = " << chi2_value << std::endl;
    }

    // 关闭文件
    file->Close();
    return 0;
}