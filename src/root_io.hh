#ifndef root_io_h
#define root_io_h 1

#include <iostream>
#include <fstream>
#include <any>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>
#include <filesystem>
namespace fs = std::filesystem;

#include "TFile.h"
#include "TString.h"
#include "TChain.h"
#include "TTree.h"

struct branch
{
    std::string name, type;
    std::string leaflist = "";
    int index = -1;
    void *value;
    void autofill();
};

struct tree
{
    TTree *ttree;
    std::vector<std::string> branch_names;
    std::map<std::string, branch> branches;
};

class RootReader
{
protected:
    static const int MAX_MULTI = 128;
    std::vector<int> addr_int;
    std::vector<double> addr_double;
    std::vector<std::array<int, MAX_MULTI>> addr_aint;
    std::vector<std::array<double, MAX_MULTI>> addr_adouble;

public:
    RootReader(const std::string &tr_name);
    RootReader(const std::string &path, const std::string &tr_name);
    ~RootReader();

    std::vector<fs::path> filenames;

    void add_file(const std::string &path);
    void set_branches(const std::vector<branch> &branches);

    TChain *tree;
    std::vector<std::string> branch_names;
    std::unordered_map<std::string, branch> branches;
    std::map<std::string, std::any> get_entry(int iEvt);
};

class RootWriter
{
protected:
    static const int MAX_MULTI = 128;
    std::vector<int> addr_int;
    std::vector<double> addr_double;
    std::vector<std::array<int, MAX_MULTI>> addr_aint;
    std::vector<std::array<double, MAX_MULTI>> addr_adouble;

public:
    RootWriter(const std::string &path, const std::string &tr_name, const std::string &option = "RECREATE");
    RootWriter(const std::string &path, const std::initializer_list<std::string> &tr_names, const std::string &option = "RECREATE");
    ~RootWriter();

    fs::path path;
    TFile *file;
    std::map<std::string, tree> trees;
    void set_branches(const std::string &tr_name, const std::vector<branch> &branches);

    void set(const std::string &tr_name, const std::string &br_name, const void *source, std::size_t nbytes);
    void set(const std::string &tr_name, const std::string &br_name, int source);
    void set(const std::string &tr_name, const std::string &br_name, double source);
    void set(const std::string &tr_name, const std::string &br_name, std::vector<int> &source);
    void set(const std::string &tr_name, const std::string &br_name, std::vector<double> &source);
    void set(const std::string &tr_name, const std::string &br_name, int size, int *source);
    void set(const std::string &tr_name, const std::string &br_name, int size, double *source);

    void fill();
    void fill(const std::string &tr_name);
    void write();
};

#endif
