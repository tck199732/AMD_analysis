#include "root_io.hh"

void branch::autofill()
{
    if (this->type.back() == ']')
    {
        int pos = std::find(this->type.begin(), this->type.end(), '[') - this->type.begin();
        std::string size_char = this->type.substr(pos + 1, this->type.length() - pos - 2);

        if (size_char == "")
        {
            throw std::invalid_argument("empty array size is not allowed in leaflist.");
        }

        char type_char = toupper(this->type[0]);
        this->leaflist = Form("%s[%s]/%c", this->name.c_str(), size_char.c_str(), type_char);
        this->type = this->type.substr(0, pos) + "[]";
    }
    else
    {
        char type_char = toupper(this->type[0]);
        this->leaflist = Form("%s/%c", this->name.c_str(), type_char);
    }
    return;
}

RootReader::RootReader(const std::string &tr_name)
{
    this->tree = new TChain(tr_name.c_str());
}

RootReader::RootReader(const std::string &path, const std::string &tr_name)
{
    this->tree = new TChain(tr_name.c_str());

    fs::path pth = path;
    if (fs ::exists(pth))
    {
        this->filenames.push_back(path);
        this->add_file(path);
    }
}

RootReader::~RootReader()
{
    delete this->tree;
}

void RootReader::add_file(const std::string &path)
{

    this->tree->Add(path.c_str());
}

void RootReader::set_branches(const std::vector<branch> &branches)
{
    for (auto &branch : branches)
    {
        this->branch_names.push_back(branch.name);
        this->branches[branch.name] = branch;
    }

    this->addr_int.resize(this->branches.size());
    this->addr_double.resize(this->branches.size());
    this->addr_aint.resize(this->branches.size());
    this->addr_adouble.resize(this->branches.size());
    this->tree->SetBranchStatus("*", false);

    int index = 0;
    for (auto &br_name : this->branch_names)
    {
        branch *branch = &this->branches[br_name];
        branch->index = index;

        this->tree->SetBranchStatus(br_name.c_str(), true);

        if (branch->type == "int")
        {
            this->tree->SetBranchAddress(br_name.c_str(), &this->addr_int[branch->index]);
            branch->value = &this->addr_int[branch->index];
        }
        else if (branch->type == "double")
        {
            this->tree->SetBranchAddress(br_name.c_str(), &this->addr_double[branch->index]);
            branch->value = &this->addr_double[branch->index];
        }
        else if (branch->type == "int[]")
        {
            this->tree->SetBranchAddress(br_name.c_str(), &this->addr_aint[branch->index][0]);
            branch->value = &this->addr_aint[branch->index][0];
        }
        else if (branch->type == "double[]")
        {
            this->tree->SetBranchAddress(br_name.c_str(), &this->addr_adouble[branch->index][0]);
            branch->value = &this->addr_adouble[branch->index][0];
        }

        index++;
    }
    return;
}

std::map<std::string, std::any> RootReader::get_entry(int iEvt)
{

    this->tree->GetEntry(iEvt);
    std::map<std::string, std::any> buffer;
    for (auto &[name, br] : this->branches)
    {
        if (br.type == "int")
        {
            buffer[br.name] = *static_cast<int *>(br.value);
        }
        else if (br.type == "double")
        {
            buffer[br.name] = *static_cast<double *>(br.value);
        }
        else if (br.type == "int[]")
        {
            buffer[br.name] = static_cast<int *>(br.value);
        }
        else if (br.type == "double[]")
        {
            buffer[br.name] = static_cast<double *>(br.value);
        }
    }
    return buffer;
}

RootWriter::RootWriter(const std::string &path, const std::string &tr_name, const std::string &option)
{
    this->path = path;
    this->file = new TFile(path.c_str(), option.c_str());
    this->trees[tr_name].ttree = new TTree(tr_name.c_str(), tr_name.c_str());
}

RootWriter::RootWriter(const std::string &path, const std::initializer_list<std::string> &tr_names, const std::string &option)
{
    this->path = path;
    this->file = new TFile(path.c_str(), option.c_str());
    for (auto &tr_name : tr_names)
    {
        this->trees[tr_name].ttree = new TTree(tr_name.c_str(), tr_name.c_str());
    }
}

RootWriter::~RootWriter()
{
    this->file->Close();
}

void RootWriter::set_branches(const std::string &tr_name, const std::vector<branch> &branches)
{
    tree *tree = &this->trees[tr_name];
    for (auto &branch : branches)
    {
        tree->branch_names.push_back(branch.name);
        tree->branches[branch.name] = branch;
    }

    auto resize = [tree](auto &&...args)
    {
        (args.resize(tree->branches.size()), ...);
    };
    resize(this->addr_int, this->addr_double, this->addr_aint, this->addr_adouble);

    int index = 0;
    for (auto &br_name : tree->branch_names)
    {
        branch *branch = &tree->branches[br_name];
        branch->index = index;
        const char *name = branch->name.c_str();
        const char *leaflist = branch->leaflist.c_str();

        if (branch->type == "int")
        {
            tree->ttree->Branch(name, &this->addr_int[index], leaflist);
            branch->value = &this->addr_int[index];
        }
        else if (branch->type == "double")
        {
            tree->ttree->Branch(name, &this->addr_double[index], leaflist);
            branch->value = &this->addr_double[index];
        }

        else if (branch->type == "int[]")
        {
            tree->ttree->Branch(name, &this->addr_aint[index][0], leaflist);
            branch->value = &this->addr_aint[index][0];
        }

        else if (branch->type == "double[]")
        {
            tree->ttree->Branch(name, &this->addr_adouble[index][0], leaflist);
            branch->value = &this->addr_adouble[index][0];
        }

        index++;
    }
    return;
}

void RootWriter::set(const std::string &tr_name, const std::string &br_name, const void *source, std::size_t nbytes)
{
    std::memcpy(this->trees[tr_name].branches[br_name].value, source, nbytes);
}

void RootWriter::set(const std::string &tr_name, const std::string &br_name, int source)
{
    this->set(tr_name, br_name, &source, sizeof(source));
}

void RootWriter::set(const std::string &tr_name, const std::string &br_name, double source)
{
    this->set(tr_name, br_name, &source, sizeof(source));
}

void RootWriter::set(const std::string &tr_name, const std::string &br_name, std::vector<int> &source)
{
    this->set(tr_name, br_name, &source[0], sizeof(source[0]) * source.size());
}

void RootWriter::set(const std::string &tr_name, const std::string &br_name, std::vector<double> &source)
{
    this->set(tr_name, br_name, &source[0], sizeof(source[0]) * source.size());
}

void RootWriter::set(const std::string &tr_name, const std::string &br_name, int size, int *source)
{
    this->set(tr_name, br_name, source, sizeof(source[0]) * size);
}

void RootWriter::set(const std::string &tr_name, const std::string &br_name, int size, double *source)
{
    this->set(tr_name, br_name, source, sizeof(source[0]) * size);
}

void RootWriter::fill()
{
    for (auto &[tr_name, tr] : this->trees)
    {
        this->trees[tr_name].ttree->Fill();
    }
}

void RootWriter::fill(const std::string &tr_name)
{
    this->trees[tr_name].ttree->Fill();
}

void RootWriter::write()
{
    this->file->cd();
    this->file->Write();
}
