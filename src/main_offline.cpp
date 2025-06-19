// g++ -std=c++17 -O3 -march=native -fopenmp -I../include -I/usr/include/eigen3 -I/home/kim/EigenRand offlineMppiEstimation_2tag.cpp main_offline.cpp -o mppi
// ./mppi
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "offlineMppiEstimation_2tag.h"

bool loadImu(const std::string& path, std::vector<ImuData>& imuAll) {
    std::ifstream ifs(path);
    if(!ifs) { std::cerr << "Cannot open " << path << std::endl; return false; }
    std::string line;
    while(std::getline(ifs, line)) {
        if(line.empty()) continue;
        std::istringstream iss(line);
        ImuData d;
        iss >> d.timeStamp
            >> d.acc.x() >> d.acc.y() >> d.acc.z()
            >> d.gyr.x() >> d.gyr.y() >> d.gyr.z();
        imuAll.push_back(d);
    }
    return true;
}

bool loadUwb(const std::string& path,
             std::vector<UwbData>& uwbLeft,
             std::vector<UwbData>& uwbRight)
{
    std::ifstream ifs(path);
    if(!ifs) { std::cerr << "Cannot open " << path << std::endl; return false; }
    std::string line;
    while(std::getline(ifs, line)) {
        if(line.empty()) continue;
        std::istringstream iss(line);
        UwbData d;
        iss >> d.id >> d.timeStamp;
        d.ranges.resize(8);
        for(int i = 0; i < 8; ++i) iss >> d.ranges(i);
        if(d.id == 0) uwbLeft.push_back(d);
        else          uwbRight.push_back(d);
    }
    return true;
}

bool getNextBatch(int T,
    const std::vector<ImuData>& imuAll, size_t& imuIdx,
    const std::vector<UwbData>& uwbLeft, size_t& ulIdx,
    const std::vector<UwbData>& uwbRight, size_t& urIdx,
    std::vector<ImuData>& imuBatch,
    std::vector<UwbData>& uwbBatch)
{
    imuBatch.clear();
    uwbBatch.clear();
    while((int)uwbBatch.size() < T) {
        if(ulIdx >= uwbLeft.size() || urIdx >= uwbRight.size()) return false;
        double tsL = uwbLeft[ulIdx].timeStamp;
        double tsR = uwbRight[urIdx].timeStamp;
        if(std::abs(tsL - tsR) > 1e-6) {
            if(tsL < tsR) ++ulIdx;
            else           ++urIdx;
            continue;
        }
        double ts = tsL;
        while(imuIdx + 1 < imuAll.size() && imuAll[imuIdx+1].timeStamp < ts) ++imuIdx;
        if(imuIdx + 1 >= imuAll.size()) return false;
        const ImuData &a = imuAll[imuIdx];
        const ImuData &b = imuAll[imuIdx+1];
        double w = (ts - a.timeStamp) / (b.timeStamp - a.timeStamp);
        ImuData di;
        di.timeStamp = ts;
        di.acc = a.acc * (1-w) + b.acc * w;
        di.gyr = a.gyr * (1-w) + b.gyr * w;
        imuBatch.push_back(di);
        UwbData du;
        du.timeStamp = ts;
        du.id = 0;
        du.ranges.resize(16);
        du.ranges.head(8) = uwbLeft[ulIdx].ranges;
        du.ranges.tail(8) = uwbRight[urIdx].ranges;
        uwbBatch.push_back(du);
        ++ulIdx;
        ++urIdx;
    }
    return true;
}

int main() {
    const std::string imuPath     = "../data/imu.txt";
    const std::string uwbPath    = "../data/uwb.txt";
    const std::string anchorsPath= "../config/anchors.txt";

    std::vector<ImuData> imuAll;
    std::vector<UwbData> uwbLeft, uwbRight;
    if(!loadImu(imuPath, imuAll) || !loadUwb(uwbPath, uwbLeft, uwbRight))
        return -1;

    double imuOffset = imuAll.front().timeStamp;
    for(auto &d : imuAll) d.timeStamp -= imuOffset;
    double uwbOffset = std::min(uwbLeft.front().timeStamp, uwbRight.front().timeStamp);
    for(auto &d : uwbLeft)  d.timeStamp -= uwbOffset;
    for(auto &d : uwbRight) d.timeStamp -= uwbOffset;

    Eigen::Matrix<double,3,8> anchors;
    std::ifstream ifs(anchorsPath);
    for(int i = 0; i < 8; ++i)
        ifs >> anchors(0,i) >> anchors(1,i) >> anchors(2,i);

    mppiEstimation mppi;
    mppi.setAnchorPositions(anchors);

    std::ofstream poseFile("../data/mppi_pose.txt");
    if(!poseFile) return -1;
    poseFile << "# time x y z\n";

    size_t imuIdx = 0, ul = 0, ur = 0;
    double prevTs = 0;

    while(true) {
        std::vector<ImuData> imuBatch;
        std::vector<UwbData> uwbBatch;
        if(!getNextBatch(mppi.T, imuAll, imuIdx,
                         uwbLeft, ul, uwbRight, ur,
                         imuBatch, uwbBatch)) break;

        double ts = uwbBatch.front().timeStamp;
        if(prevTs > 0) mppi.setDt(ts - prevTs);
        prevTs = ts;

        mppi.solve(uwbBatch, imuBatch);
        mppi.move(imuBatch.back());
        STATE st = mppi.xInit;
        poseFile << ts << " "
        << st.p.x() << " "
        << st.p.y() << " "
        << st.p.z() << "\n";
        std::cout<< ts <<" pos=["<< st.p.transpose() <<"]\n";
    }
    poseFile.close();

    int ret = std::system("python3 ../scripts/plot.py"); (void)ret;

    return 0;
}
