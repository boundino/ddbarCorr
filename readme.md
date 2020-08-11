ddbarCorr
====

* Get the repository
```
git clone https://github.com/boundino/ddbarCorr.git --recursive
```

* Run the macro
```
cd ddbarCorr/mainana/
./run_ddbar.sh 1 1
```

* Run 2D fit using cmake
```
cd ddbarCorr/mainana/
./tree.sh 1 1
```
Note that cmake doesn't work properly with `/raid5/root/root-v6.16.00-gcc731/`.
This is becasue the Cmake files of that installation isn't set up properly:
```grep Macros /raid5/root/root-v6.16.00-gcc731/ROOTUseFile.cmake``` points to a file with no read permission.
To work around that, switch to another installation. E.g. ```locate thisroot.sh``` and source either the root 6.20 or 6.22 in the output.
