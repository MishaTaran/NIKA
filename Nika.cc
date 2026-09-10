#include <Loader.hh>

int main(int argc, char **argv) {

    //G4cout.rdbuf(nullptr);
    Loader *loader = new Loader(argc, argv);
    delete loader;
    return 0;
}
