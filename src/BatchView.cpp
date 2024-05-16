#include "Batchview.h"

BatchView::BatchView(sf::RenderWindow& window, Landscape* pLand, Community* pCommunity) : 
	pLandscape{pLand}, 
	pComm{pCommunity} {
	
	// Set cell size such that both dimensions fit on screen
	dimX = pLandscape->getLandParams().dimX;
	dimY = pLandscape->getLandParams().dimY;
	cellSize = min(1920u / dimX, 1080u / dimY);

	window.create(sf::VideoMode{dimX * cellSize, dimY * cellSize}, "RangeShifter Batch");
	window.setFramerateLimit(144);
}


// ---------------------------------------
// Collect and process user window input
// e.g. clicking, scrolling, typing etc.
// ---------------------------------------
void BatchView::collectUserInput(sf::RenderWindow& window) {
	if (window.isOpen()) {
		for (auto event = sf::Event{}; window.pollEvent(event);)
		{
			// Close the window
			if (event.type == sf::Event::Closed)
			{
				window.close();
			}
			// else, add more events...
		}
	}
}

void BatchView::drawLandscape(sf::RenderWindow& window) {

	for (int x = 0; x < dimX; x++) {
		for (int y = 0; y < dimY; y++) {
			Cell* pCell = pLandscape->findCell(x, y);
			if (pCell != 0) {
				sf::RectangleShape c(sf::Vector2f(cellSize, cellSize));
				c.setPosition(cellSize * x, cellSize * y);
				c.setFillColor(sf::Color(pRandom->IRandom(0, 255), pRandom->IRandom(0, 255), pRandom->IRandom(0, 255)));
				window.draw(c);
			}
		}
	}

}

void BatchView::drawCommunity(sf::RenderWindow& window) {

	// Erase previous community
	window.clear();
	drawLandscape(window);

	// Draw current community
	// ...

	window.display();
}