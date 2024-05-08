#include "Batchview.h"

BatchView::BatchView(Landscape* pLand, Community* pCommunity) : pLandscape{pLand}, pComm{pCommunity} {

	// Open a window
	auto window = sf::RenderWindow{ { 1920u, 1080u }, "RangeShifter Batch" };
	window.setFramerateLimit(144);
	window.display();
	/*
	while (window.isOpen())
	{
		window.clear();
		window.display();
	} */
}


// ---------------------------------------
// Collect and process user window input
// e.g. clicking, scrolling, typing etc.
// ---------------------------------------
void BatchView::collectUserInput() {
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

void BatchView::drawLandscape() {

}

void BatchView::drawCommunity() {
	window.clear();
	window.display();
}