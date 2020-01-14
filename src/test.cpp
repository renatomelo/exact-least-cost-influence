
#include <iostream>
#include <vector>

class Phi
{
    static int _index;

    private:
        std::vector<int> nodes;
        int index;
    
    public:
        Phi() {
            index = _index++;
        }

        operator int() {
            return index;
        }
};

int Phi::_index = 0;

int main(){

    int v[4] = {3, 6, 9, 12};
    for (int i = 0; i < 4; i++) {
        Phi phi;
        std::cout << v[phi] << std::endl;
    }

}