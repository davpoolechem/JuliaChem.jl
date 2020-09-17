#/bin/bash

#== download libint code generator ==#
cd /home/travis 
git clone --quiet https://github.com/evaleev/libint.git

#== patch libint Makefiles to remove CMake references ==#
cd libint
echo "33c33" &> Makefile.patch
echo "< install:: all install_pkgconfig install_cmake install_inc install_data" >> Makefile.patch
echo "---" >> Makefile.patch
echo "> install:: all install_pkgconfig install_inc install_data" >> Makefile.patch
echo "43,48d42" >> Makefile.patch
echo "< endif" >> Makefile.patch
echo "< " >> Makefile.patch
echo "< ifdef cmakedir" >> Makefile.patch
echo "< install_cmake::" >> Makefile.patch
echo -e "< \t\$(INSTALL) \$(INSTALLDIROPT) \$(DESTDIR)\$(cmakedir)" >> Makefile.patch
echo -e "< \t\$(INSTALL) \$(INSTALLLIBOPT) \$(SRCTOPDIR)/FindLibint2.cmake \$(DESTDIR)\$(cmakedir)" >> Makefile.patch
patch Makefile Makefile.patch

cd src/lib
echo "2c2" &> MakeRules.in.patch
echo "< .PHONY: default export install install_inc install_pkgconfig install_cmake install_target clean oclean distclean targetclean realclean" >> MakeRules.in.patch
echo "---" >> MakeRules.in.patch
echo "> .PHONY: default export install install_inc install_pkgconfig install_target clean oclean distclean targetclean realclean" >> MakeRules.in.patch
echo "21c21" >> MakeRules.in.patch
echo "< install:: install_inc install_target install_pkgconfig install_cmake install_data" >> MakeRules.in.patch
echo "---" >> MakeRules.in.patch
echo "> install:: install_inc install_target install_pkgconfig install_data" >> MakeRules.in.patch
echo "51,56d50" >> MakeRules.in.patch
echo "< endif" >> MakeRules.in.patch
echo "< " >> MakeRules.in.patch
echo "< ifdef cmakedir" >> MakeRules.in.patch
echo "< install_cmake::" >> MakeRules.in.patch
echo -e "< \t\$(INSTALL) \$(INSTALLDIROPT) \$(DESTDIR)\$(cmakedir)" >> MakeRules.in.patch
echo -e "< \t\$(INSTALL) \$(INSTALLLIBOPT) \$(SRCTOPDIR)/FindLibint2.cmake \$(DESTDIR)\$(cmakedir)" >> MakeRules.in.patch
patch MakeRules.in MakeRules.in.patch

#== install libint ==#
cd /home/travis/libint
./configure \
  --enable-shared=yes --prefix=/home/travis/libint-install \
  --enable-1body=0 --enable-eri=0 --with-max-am=3 \
  --with-multipole-max-order=3

make -j2 -s
make -j2 -s install
