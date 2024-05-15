#ifndef ViewH
#define ViewH

using namespace std;

#include <SFML/Graphics.hpp>
#include "RScore/Landscape.h"
#include "RScore/Community.h"

class BatchView {
public:
	BatchView::BatchView(Landscape* pLand, Community* pCommunity);
	bool isOpen(sf::RenderWindow& window) const { return window.isOpen(); }
	//void close() { window.close(); }
	void collectUserInput(sf::RenderWindow& window);
	void drawLandscape(sf::RenderWindow& window);
	void drawCommunity(sf::RenderWindow& window);

private:
	//sf::RenderWindow window;

	Landscape* pLandscape;
	Community* pComm;
};

#endif ViewH