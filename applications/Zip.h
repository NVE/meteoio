#ifndef ZIP_H
#define ZIP_H

#include <zip.h> // https://libzip.org/documentation/libzip.html
#include <dirent.h>
#include <sys/stat.h>
#include "InternalServerError.h"
#include <iostream>

using namespace std;

// Custom request handler
class Zip {
public:
    static bool is_dir(const string& dir)
    {
        struct stat st;
        stat(dir.c_str(), &st);
        return S_ISDIR(st.st_mode);
    }

    static void zipDirectory(string &resultDir) {
        // Zip the result directory
        DIR *dir;
        int zipError = 0;
        cout << "Creating zip at " << resultDir << ".zip" << endl;
        zip_t *zipFile = zip_open((resultDir + ".zip").c_str(), ZIP_CREATE | ZIP_EXCL, &zipError);
        if(zipFile == nullptr) {
            zip_error_t ziperror;
            zip_error_init_with_code(&ziperror, zipError);
            throw InternalServerError("Error creating result zip! Code=" + string(zip_error_strerror(&ziperror)));
        }
        try {
            dir = opendir(resultDir.c_str());
            if(dir == NULL)
                throw InternalServerError("Result directory could not be opened");
            struct dirent *entry;
            while ((entry = readdir(dir))) {
                string fullname = resultDir + "/" + entry->d_name;
                if(is_dir(fullname)) continue;
                cout << "Adding file " << fullname << " to zip" << endl;
                zip_source_t *source = zip_source_file(zipFile, fullname.c_str(), 0, 0);
                if (source == nullptr) {
                    throw InternalServerError("Error reading source of file " + string(entry->d_name) + "! Code=" + string(zip_strerror(zipFile)));
                }
                if (zip_file_add(zipFile, entry->d_name, source, ZIP_FL_ENC_UTF_8) < 0) {
                    zip_source_free(source);
                    throw InternalServerError("Error adding file " + string(entry->d_name) + " to result zip! Code=" + string(zip_strerror(zipFile)));
                }
            }
        }
        catch (exception &e) {
            zip_close(zipFile);
            closedir(dir);
            throw e;
        }
        closedir(dir);
        zip_close(zipFile);
    }

};

#endif // ZIP_H