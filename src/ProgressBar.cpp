#include <iostream>
#include <string>
#include <array>

#include "TString.h"
class ProgressBar
{
public:
    int njobs;
    int njobs_done;
    std::string progress_name;
    ProgressBar(const int &njobs, const std::string &name = "Progress")
    {
        this->progress_name = name;
        this->njobs = njobs;
        this->njobs_done = 0;
    }
    void Update()
    {
        if (this->njobs_done == this->njobs)
        {
            return;
        }
        this->njobs_done++;

        if (this->njobs_done % (this->njobs / 100) == 0)
        {
            double progress = (double)njobs_done / njobs;
            int barWidth = 70;
            std::cout << Form("%14s: [", this->progress_name.c_str());
            std::cout << std::string(barWidth * progress, '=') << std::string(barWidth * (1 - progress), ' ') << "] " << int(progress * 100.0) << "%\r";
            std::cout << std::flush;
        }

        if (this->njobs_done == this->njobs)
        {
            this->Finish();
        }
    }
    void Finish()
    {
        std::cout << Form("%14s: [", this->progress_name.c_str());
        std::cout << std::string(70, '=') << "] " << 100 << "%\r";
        std::cout << std::endl;
    }
};