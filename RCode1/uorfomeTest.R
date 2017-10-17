library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("/export/valenfs/projects/uORFome/RCode1/CppFiles/main.cpp")

find_orfs_in_specific_frame("ATGATGTAATAA",
                            "ATG|TGA|GGG",
                            "TAA|AAT|ATA",
                            T,
                            0)