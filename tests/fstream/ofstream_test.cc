#include <stdlib.h>
#include <meteoio/MeteoIO.h>
#include <sys/stat.h>
#include <fstream>
#include <meteoio/FStream.h>


using namespace mio; //The MeteoIO namespace is called mio

int main() {
    std::cerr<< "Testing the generation of directories from path"<<std::endl;
    std::string filename1 = "./this/is/a/test.txt";
    std::string filename2 = "./this/is/a/second_test.txt";
    ofilestream ofs1(filename1.c_str());
    ofilestream ofs2(filename2.c_str());
    std::string instr ="This is a test";
    ofs1 << instr;
    ofs2 << instr;
    ofs1.close();
    ofs2.close();
	struct stat sb;

    bool dir_status = false;
    if (stat(mio::FileUtils::getPath(filename1,false).c_str(), &sb) == 0) {
        dir_status = true;
    } else {
        throw IOException("Output directory was not created.");
    }

    std::ifstream infile1(filename1);
    std::ifstream infile2(filename2);
    std::string item;
    bool filestatus = false;
    while (getline(infile1, item)) {
        if (item != instr) {
            throw IOException("File was not created, or written to!");
        } else {
            filestatus = true;
        }
    }
    while (getline(infile2, item)) {
        if (item != instr) {
            throw IOException("File was not created, or written to!");
        } else {
            filestatus = true;
        }
    }
    std::cerr << "Everything works as expected"<< std::endl;
    std::cerr << "----------------------------" << std::endl;

#ifdef LIMIT_WRITE_ACCESS
    std::cerr << "Testing to write a not allowed directory with limited writing access"<<std::endl;
    std::string filename3 ="../../this_shouldnt_be_outside_of_fstream/works.txt";
    ofilestream file3(filename3);
    file3 << "hiii";
    file3.close();
    if (FileUtils::directoryExists("../../this_shouldnt_be_outside_of_fstream")) {
        std::cerr << "limiting write acces did not work"<<std::endl;
    } else if (FileUtils::directoryExists("this_shouldnt_be_outside_of_fstream")) {
        std::cerr << "limiting write acess works perfectly"<< std::endl;
    } else {
        std::cerr << "directory was not created at all in write acces limitation scope"<< std::endl;
    }
    int a = system ("rm -rf this_shouldnt_be_outside_of_fstream");
    std::cerr << "----------------------------" << std::endl;
#endif

    int _ = system("rm -rf this");
    std::cerr << "Checking for calls to std::ofstream outside of wrapper"<< std::endl;
    int sys_o = system("rgrep --exclude=\"*.cc.o\" \"ofstream\" ../../meteoio | grep -vE \"^../../meteoio/FStream.\" | grep -vE \"ofstream::app\" |grep -vE \"ofstream::out\">tmp_ofstream_instances.txt");
    std::ifstream ofstream_instances("tmp_ofstream_instances.txt");
    while (getline(ofstream_instances, item,':')) {
        if (item.empty()) continue;
        else {
            std::cerr <<"Instance of std::ofstream found in file " <<item;
            exit(1);
        }
    }
    _ = system("rm tmp_ofstream_instances.txt");
    std::cerr << "Nothing found"<< std::endl;

    return !filestatus || !dir_status;
}