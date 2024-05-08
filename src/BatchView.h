#ifndef ViewH
#define ViewH

using namespace std;

#include <SFML/Graphics.hpp>
#include "RScore/Landscape.h"
#include "RScore/Community.h"

class BatchView {
public:
	BatchView::BatchView(Landscape* pLand, Community* pCommunity);
	bool isOpen() const { return window.isOpen(); }
	void close() { window.close(); }
	void collectUserInput();
	void drawLandscape();
	void drawCommunity();

private:
	sf::RenderWindow window;

	Landscape* pLandscape;
	Community* pComm;
};

#endif ViewH