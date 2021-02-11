if [ ! -f JuliaChem.o ]; then
  echo "Creating JuliaChem object file..."
  julia --startup-file=no --output-o JuliaChem.o --optimize=3 --check-bounds=no --math-mode=fast --inline=yes --compiled-modules=yes -J"$JULIA_ROOT/lib/julia/sys.so" custom_sysimg.jl
fi

echo "Creating JuliaChem system image..." 
gcc -Ofast -shared -o JuliaChem.so -Wl,--whole-archive JuliaChem.o -Wl,--no-whole-archive -L"$JULIA_ROOT/lib" -ljulia
