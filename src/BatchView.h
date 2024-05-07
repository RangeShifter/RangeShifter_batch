#ifndef ViewH
#define ViewH

using namespace std;

#include <SFML/Graphics.hpp>
#include "RScore/Landscape.h"
#include "RScore/Community.h"

class BatchView {
public:
	BatchView::BatchView();

private:
	sf::RenderWindow window;

	Landscape* pLandscape;
	Community* pComm;
};

#endif ViewH