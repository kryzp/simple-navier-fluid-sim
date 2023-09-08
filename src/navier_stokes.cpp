#include <iostream>
#include <SFML/Graphics.hpp>
#include <cmath>
#include <vector>

template <typename T>
inline void util_swap(T& x, T& y)
{
	T tmp = x;
	x = y;
	y = tmp;
}

#define INDEX(x, y) (((y)*MAP_WIDTH)+(x))
#define MIN(x, y) ((x > y) ? y : x)
#define MAX(x, y) ((x > y) ? x : y)
#define WINDOW_WIDTH 720
#define WINDOW_HEIGHT 720
#define SCALE 8
#define MAP_WIDTH (WINDOW_WIDTH / SCALE)
#define MAP_HEIGHT (WINDOW_HEIGHT / SCALE)
#define SWAP(_x, _y) (util_swap((_x), (_y)))
#define ITERATIONS 50
#define VISCOSITY 0.0
#define DIFFUSION_RATE 0.0

struct Cell
{
	sf::Vector2<double> velocity;
	double density;
};

Cell* cells = nullptr;
Cell* cells_prev = nullptr;
double* divergence = nullptr;
double* p = nullptr;

double calculate_curl(int x, int y)
{
	int left = INDEX(x - 1, y);
	int right = INDEX(x + 1, y);
	int top = INDEX(x, y - 1);
	int bottom = INDEX(x, y + 1);

	double dxdy = (cells[bottom].velocity.x - cells[top].velocity.x);
	double dydx = (cells[left].velocity.y - cells[right].velocity.y);

	return dxdy - dydx;
}

void boundary_conditions()
{
	//cells[INDEX(25, 25)].velocity.x = 0.0;
	//cells[INDEX(25, 25)].velocity.y = 1.0;
	//cells[INDEX(25, 25)].density = 100.0;

	for (int i = 0; i < MAP_WIDTH; i++)
	{
		cells[INDEX(i, 0)].velocity.y = -cells[INDEX(i, 1)].velocity.y;
		cells[INDEX(i, MAP_HEIGHT - 1)].velocity.y = -cells[INDEX(i, MAP_HEIGHT - 2)].velocity.y;

		cells[INDEX(i, 0)].density = cells[INDEX(i, 1)].density;
		cells[INDEX(i, MAP_HEIGHT - 1)].density = cells[INDEX(i, MAP_HEIGHT - 2)].density;

		divergence[INDEX(i, 0)] = divergence[INDEX(i, 1)];
		divergence[INDEX(i, MAP_HEIGHT - 1)] = divergence[INDEX(i, MAP_HEIGHT - 2)];

		p[INDEX(i, 0)] = p[INDEX(i, 1)];
		p[INDEX(i, MAP_HEIGHT - 1)] = p[INDEX(i, MAP_HEIGHT - 2)];
	}

	for (int i = 0; i < MAP_HEIGHT; i++)
	{
		cells[INDEX(0, i)].velocity.x = -cells[INDEX(1, i)].velocity.x;
		cells[INDEX(MAP_WIDTH - 1, i)].velocity.x = -cells[INDEX(MAP_WIDTH - 2, i)].velocity.x;

		cells[INDEX(0, i)].density = cells[INDEX(1, i)].density;
		cells[INDEX(MAP_WIDTH - 1, i)].density = cells[INDEX(MAP_WIDTH - 2, i)].density;

		divergence[INDEX(0, i)] = divergence[INDEX(1, i)];
		divergence[INDEX(MAP_WIDTH - 1, i)] = divergence[INDEX(MAP_WIDTH - 2, i)];

		p[INDEX(0, i)] = p[INDEX(1, i)];
		p[INDEX(MAP_WIDTH - 1, i)] = p[INDEX(MAP_WIDTH - 2, i)];
	}
}

void diffuse_density(double dt)
{
	double a = dt * DIFFUSION_RATE;

	for (int k = 0; k < ITERATIONS; k++)
	{
		for (int y = 1; y < MAP_HEIGHT - 1; y++)
		{
			for (int x = 1; x < MAP_WIDTH - 1; x++)
			{
				int i = INDEX(x, y);
				int top = INDEX(x, y - 1);
				int bottom = INDEX(x, y + 1);
				int left = INDEX(x - 1, y);
				int right = INDEX(x + 1, y);

				cells[i].density = cells_prev[i].density;

				/*
				cells[i].density = (cells_prev[i].density + a*(
					cells[left].density +
					cells[right].density +
					cells[top].density +
					cells[bottom].density
				)) / (1.0 + 4.0*a);
				 */
			}
		}

		boundary_conditions();
	}
}

void advect_density(double dt)
{
	for (int i = 1; i < MAP_WIDTH - 1; i++)
	{
		for (int j = 1; j < MAP_HEIGHT - 1; j++)
		{
			double x = i - dt * (double)MAP_WIDTH  * cells_prev[INDEX(i, j)].velocity.x;
			double y = j - dt * (double)MAP_HEIGHT * cells_prev[INDEX(i, j)].velocity.y;

			if (x < 0.5)
				x = 0.5;

			if (y < 0.5)
				y = 0.5;

			if (x > (double)MAP_WIDTH - 0.5)
				x = (double)MAP_WIDTH - 0.5;

			if (y > (double)MAP_HEIGHT - 0.5)
				y = (double)MAP_HEIGHT - 0.5;

			int i0 = (int)x;
			int j0 = (int)y;
			int i1 = i0 + 1;
			int j1 = j0 + 1;

			double s1 = x - i0;
			double s0 = 1 - s1;
			double t1 = y - j0;
			double t0 = 1 - t1;

			cells[INDEX(i, j)].density =
				s0 * (t0 * cells_prev[INDEX(i0, j0)].density + t1 * cells_prev[INDEX(i0, j1)].density) +
				s1 * (t0 * cells_prev[INDEX(i1, j0)].density + t1 * cells_prev[INDEX(i1, j1)].density);
		}
	}

	boundary_conditions();
}

void density_tick(double dt)
{
	SWAP(cells_prev, cells);
	diffuse_density(dt);

	SWAP(cells_prev, cells);
	advect_density(dt);
}

void project_velocity()
{
	for (int i = 1; i < MAP_WIDTH - 1; i++)
	{
		for (int j = 1; j < MAP_HEIGHT - 1; j++)
		{
			divergence[INDEX(i, j)] =
				-0.5 / (double)MAP_WIDTH  * (cells[INDEX(i + 1, j)].velocity.x - cells[INDEX(i - 1, j)].velocity.x) +
				-0.5 / (double)MAP_HEIGHT * (cells[INDEX(i, j + 1)].velocity.y - cells[INDEX(i, j - 1)].velocity.y);

			p[INDEX(i, j)] = 0.0;
		}
	}

	boundary_conditions();

	for (int k = 0; k < ITERATIONS; k++)
	{
		for (int i = 1; i < MAP_WIDTH - 1; i++)
		{
			for (int j = 1; j < MAP_HEIGHT - 1; j++)
			{
				p[INDEX(i, j)] = (divergence[INDEX(i, j)] +
					p[INDEX(i - 1, j)] + p[INDEX(i + 1, j)] +
					p[INDEX(i, j - 1)] + p[INDEX(i, j + 1)]) / 4.0;
			}
		}

		boundary_conditions();
	}

	for (int i = 1; i < MAP_WIDTH - 1; i++)
	{
		for (int j = 1; j < MAP_HEIGHT - 1; j++)
		{
			cells[INDEX(i, j)].velocity.x -= 0.5 * (p[INDEX(i + 1, j)] - p[INDEX(i - 1, j)]) * (double)MAP_WIDTH;
			cells[INDEX(i, j)].velocity.y -= 0.5 * (p[INDEX(i, j + 1)] - p[INDEX(i, j - 1)]) * (double)MAP_HEIGHT;
		}
	}

	boundary_conditions();
}

void diffuse_velocity(double dt)
{
	double a = dt * VISCOSITY;

	for (int k = 0; k < ITERATIONS; k++)
	{
		for (int y = 1; y < MAP_HEIGHT - 1; y++)
		{
			for (int x = 1; x < MAP_WIDTH - 1; x++)
			{
				int i = INDEX(x, y);
				int top = INDEX(x, y - 1);
				int bottom = INDEX(x, y + 1);
				int left = INDEX(x - 1, y);
				int right = INDEX(x + 1, y);

				cells[i].velocity = (cells_prev[i].velocity + a*(
					cells[left].velocity / 4.0 +
					cells[right].velocity / 4.0 +
					cells[top].velocity / 4.0 +
					cells[bottom].velocity / 4.0
				)) / (1.0 + a);
			}
		}

		boundary_conditions();
	}
}

void advect_velocity(double dt)
{
	for (int i = 1; i < MAP_WIDTH - 1; i++)
	{
		for (int j = 1; j < MAP_HEIGHT - 1; j++)
		{
			double x = i - dt * (double)MAP_WIDTH  * cells_prev[INDEX(i, j)].velocity.x;
			double y = j - dt * (double)MAP_HEIGHT * cells_prev[INDEX(i, j)].velocity.y;

			if (x < 0.5)
				x = 0.5;

			if (y < 0.5)
				y = 0.5;

			if (x > (double)MAP_WIDTH - 0.5)
				x = (double)MAP_WIDTH - 0.5;

			if (y > (double)MAP_HEIGHT - 0.5)
				y = (double)MAP_HEIGHT - 0.5;

			int i0 = (int)x;
			int j0 = (int)y;
			int i1 = i0 + 1;
			int j1 = j0 + 1;

			double s1 = x - i0;
			double s0 = 1 - s1;
			double t1 = y - j0;
			double t0 = 1 - t1;

			cells[INDEX(i, j)].velocity =
				s0 * (t0 * cells_prev[INDEX(i0, j0)].velocity + t1 * cells_prev[INDEX(i0, j1)].velocity) +
				s1 * (t0 * cells_prev[INDEX(i1, j0)].velocity + t1 * cells_prev[INDEX(i1, j1)].velocity);
		}
	}

	boundary_conditions();
}

void velocity_tick(double dt)
{
	SWAP(cells_prev, cells);
	diffuse_velocity(dt);

	project_velocity();

	SWAP(cells_prev, cells);
	advect_velocity(dt);

	project_velocity();
}

void update(double dt)
{
	velocity_tick(dt);
	density_tick(dt);

	memcpy(cells_prev, cells, MAP_WIDTH * MAP_HEIGHT * sizeof(Cell));
}

struct GradientSettings
{
	double ar, ag, ab;
	double br, bg, bb;
	double cr, cg, cb;
	double dr, dg, db;
};

sf::Color gradient(float t, GradientSettings s)
{
	double mult = pow(t, 0.1);

	double colr = mult * (s.ar + (s.br * cos(6.28318 * ((s.cr * t) + s.dr))));
	double colg = mult * (s.ag + (s.bg * cos(6.28318 * ((s.cg * t) + s.dg))));
	double colb = mult * (s.ab + (s.bb * cos(6.28318 * ((s.cb * t) + s.db))));

	colr = MIN(1.0, MAX(0.0, colr));
	colg = MIN(1.0, MAX(0.0, colg));
	colb = MIN(1.0, MAX(0.0, colb));

	return sf::Color(
		(uint8_t)(colr * 255.0),
		(uint8_t)(colg * 255.0),
		(uint8_t)(colb * 255.0)
	);
}

void smoothify_blocks(sf::Vertex* pixels, const sf::Vertex* prev_pixels)
{
	for (int y = 1; y < MAP_HEIGHT - 1; y++)
	{
		for (int x = 1; x < MAP_WIDTH - 1; x++)
		{
			pixels[4*INDEX(x, y) + 0].color.r = (
				(int)prev_pixels[4*INDEX(x, y) + 0].color.r +
				(int)prev_pixels[4*INDEX(x - 1, y) + 1].color.r +
				(int)prev_pixels[4*INDEX(x, y - 1) + 3].color.r +
				(int)prev_pixels[4*INDEX(x - 1, y - 1) + 2].color.r
			) / 4;

			pixels[4*INDEX(x, y) + 1].color.r = (
				(int)prev_pixels[4*INDEX(x, y) + 1].color.r +
				(int)prev_pixels[4*INDEX(x + 1, y) + 0].color.r +
				(int)prev_pixels[4*INDEX(x + 1, y - 1) + 3].color.r +
				(int)prev_pixels[4*INDEX(x, y - 1) + 2].color.r
			) / 4;

			pixels[4*INDEX(x, y) + 2].color.r = (
				(int)prev_pixels[4*INDEX(x, y) + 2].color.r +
				(int)prev_pixels[4*INDEX(x + 1, y) + 3].color.r +
				(int)prev_pixels[4*INDEX(x + 1, y + 1) + 0].color.r +
				(int)prev_pixels[4*INDEX(x, y + 1) + 1].color.r
			) / 4;

			pixels[4*INDEX(x, y) + 3].color.r = (
				(int)prev_pixels[4*INDEX(x, y) + 3].color.r +
				(int)prev_pixels[4*INDEX(x - 1, y) + 2].color.r +
				(int)prev_pixels[4*INDEX(x - 1, y + 1) + 1].color.r +
				(int)prev_pixels[4*INDEX(x, y + 1) + 0].color.r
			) / 4;

			pixels[4*INDEX(x, y) + 0].color.g = (
				(int)prev_pixels[4*INDEX(x, y) + 0].color.g +
					(int)prev_pixels[4*INDEX(x - 1, y) + 1].color.g +
					(int)prev_pixels[4*INDEX(x, y - 1) + 3].color.g +
					(int)prev_pixels[4*INDEX(x - 1, y - 1) + 2].color.g
			) / 4;

			pixels[4*INDEX(x, y) + 1].color.g = (
				(int)prev_pixels[4*INDEX(x, y) + 1].color.g +
					(int)prev_pixels[4*INDEX(x + 1, y) + 0].color.g +
					(int)prev_pixels[4*INDEX(x + 1, y - 1) + 3].color.g +
					(int)prev_pixels[4*INDEX(x, y - 1) + 2].color.g
			) / 4;

			pixels[4*INDEX(x, y) + 2].color.g = (
				(int)prev_pixels[4*INDEX(x, y) + 2].color.g +
					(int)prev_pixels[4*INDEX(x + 1, y) + 3].color.g +
					(int)prev_pixels[4*INDEX(x + 1, y + 1) + 0].color.g +
					(int)prev_pixels[4*INDEX(x, y + 1) + 1].color.g
			) / 4;

			pixels[4*INDEX(x, y) + 3].color.g = (
				(int)prev_pixels[4*INDEX(x, y) + 3].color.g +
					(int)prev_pixels[4*INDEX(x - 1, y) + 2].color.g +
					(int)prev_pixels[4*INDEX(x - 1, y + 1) + 1].color.g +
					(int)prev_pixels[4*INDEX(x, y + 1) + 0].color.g
			) / 4;

			pixels[4*INDEX(x, y) + 0].color.b = (
				(int)prev_pixels[4*INDEX(x, y) + 0].color.b +
					(int)prev_pixels[4*INDEX(x - 1, y) + 1].color.b +
					(int)prev_pixels[4*INDEX(x, y - 1) + 3].color.b +
					(int)prev_pixels[4*INDEX(x - 1, y - 1) + 2].color.b
			) / 4;

			pixels[4*INDEX(x, y) + 1].color.b = (
				(int)prev_pixels[4*INDEX(x, y) + 1].color.b +
					(int)prev_pixels[4*INDEX(x + 1, y) + 0].color.b +
					(int)prev_pixels[4*INDEX(x + 1, y - 1) + 3].color.b +
					(int)prev_pixels[4*INDEX(x, y - 1) + 2].color.b
			) / 4;

			pixels[4*INDEX(x, y) + 2].color.b = (
				(int)prev_pixels[4*INDEX(x, y) + 2].color.b +
					(int)prev_pixels[4*INDEX(x + 1, y) + 3].color.b +
					(int)prev_pixels[4*INDEX(x + 1, y + 1) + 0].color.b +
					(int)prev_pixels[4*INDEX(x, y + 1) + 1].color.b
			) / 4;

			pixels[4*INDEX(x, y) + 3].color.b = (
				(int)prev_pixels[4*INDEX(x, y) + 3].color.b +
					(int)prev_pixels[4*INDEX(x - 1, y) + 2].color.b +
					(int)prev_pixels[4*INDEX(x - 1, y + 1) + 1].color.b +
					(int)prev_pixels[4*INDEX(x, y + 1) + 0].color.b
			) / 4;
		}
	}
}

int main()
{
	sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Navier Stokes");

	GradientSettings grad = {
		0.938,
		0.590,
		0.718,
		0.659,
		0.610,
		0.328,
		0.388,
		0.388,
		0.296,
		2.538,
		2.478,
		0.168
	};

	sf::Vertex* pixels = new sf::Vertex[4 * MAP_WIDTH * MAP_HEIGHT];
	sf::Vertex* prev_pixels = new sf::Vertex[4 * MAP_WIDTH * MAP_HEIGHT];

	for (int y = 0; y < MAP_HEIGHT; y++)
	{
		for (int x = 0; x < MAP_WIDTH; x++)
		{
			sf::Vertex tl = {}, tr = {}, bl = {}, br = {};

			float px = x * SCALE;
			float py = y * SCALE;

			tl.position = { px        , py         };
			tr.position = { px + SCALE, py         };
			bl.position = { px        , py + SCALE };
			br.position = { px + SCALE, py + SCALE };

			int idx = INDEX(x, y);

			pixels[4*idx + 0] = tl;
			pixels[4*idx + 1] = tr;
			pixels[4*idx + 2] = br;
			pixels[4*idx + 3] = bl;
		}
	}

	cells = new Cell[MAP_WIDTH * MAP_HEIGHT];
	cells_prev = new Cell[MAP_WIDTH * MAP_HEIGHT];
	divergence = new double[MAP_WIDTH * MAP_HEIGHT];
	p = new double[MAP_WIDTH * MAP_HEIGHT];

	memset(cells, 0, MAP_WIDTH * MAP_HEIGHT * sizeof(Cell));
	memset(cells_prev, 0, MAP_WIDTH * MAP_HEIGHT * sizeof(Cell));
	memset(divergence, 0, MAP_WIDTH * MAP_HEIGHT * sizeof(double));
	memset(p, 0, MAP_WIDTH * MAP_HEIGHT * sizeof(double));

	double delta_time = 0.0;
	sf::Clock delta_clock;

	int draw_mode = 1;

	sf::Vector2i mouse_pos, prev_mouse_pos;

	while (window.isOpen())
	{
		sf::Time dt = delta_clock.restart();
		delta_time = dt.asSeconds();

		sf::Event event;
		while (window.pollEvent(event)) {
			if (event.type == sf::Event::Closed) {
				window.close();
			}
		}

		if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Left))
		{
			prev_mouse_pos = mouse_pos;
			mouse_pos = sf::Mouse::getPosition(window);

			mouse_pos /= SCALE;

			sf::Vector2i delta = mouse_pos - prev_mouse_pos;

			double len = sqrt(delta.x*delta.x + delta.y*delta.y);

			if (len >= 0.0001)
			{
				sf::Vector2<double> norm((double)delta.x / len, (double)delta.y / len);

				cells[INDEX(mouse_pos.x, mouse_pos.y)].velocity.x += norm.x * 50.0;
				cells[INDEX(mouse_pos.x, mouse_pos.y)].velocity.y += norm.y * 50.0;
			}
		}
		else
		{
			mouse_pos = prev_mouse_pos = sf::Mouse::getPosition(window);
		}

		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Space))
		{
			memset(cells, 0, MAP_WIDTH * MAP_HEIGHT * sizeof(Cell));
			memset(cells_prev, 0, MAP_WIDTH * MAP_HEIGHT * sizeof(Cell));
			memset(divergence, 0, MAP_WIDTH * MAP_HEIGHT * sizeof(double));
			memset(p, 0, MAP_WIDTH * MAP_HEIGHT * sizeof(double));

			for (int i = 0; i < MAP_WIDTH / 4; i++)
			{
				for (int j = 0; j < MAP_HEIGHT / 2; j++)
				{
					cells[INDEX(MAP_WIDTH/4 + i, MAP_HEIGHT/4 + j)].density = 0.1;
				}
			}
		}

		update(delta_time);

		for (int y = 0; y < MAP_HEIGHT; y++)
		{
			for (int x = 0; x < MAP_WIDTH; x++)
			{
				int idx = INDEX(x, y);

				if (!sf::Keyboard::isKeyPressed(sf::Keyboard::G))
				{
					double brightness = cells[idx].density * 100.0;

					//brightness = sqrt(u*u + v*v) / 50.f;
					//brightness = 1.f - (u*u + v*v) / 100.f;

					//float r = brightness * 5.f;
					//float g = brightness * 5.f;
					//float b = brightness * 5.f;

					//u = abs(u);
					//v = abs(v);

					//r += u;
					//g += v;
					//b += sqrt(u*u + v*v);

					//r = r / (1.f + r);
					//g = g / (1.f + g);
					//b = b / (1.f + b);

					//r *= 255.f;
					//g *= 255.f;
					//b *= 255.f;

					brightness = brightness / (1.0 + brightness);
					sf::Color colour = gradient(brightness, grad);

					//sf::Color colour(
					//	(uint8_t)r,
					//	(uint8_t)g,
					//	(uint8_t)b
					//);

					pixels[idx*4 + 0].color = colour;
					pixels[idx*4 + 1].color = colour;
					pixels[idx*4 + 2].color = colour;
					pixels[idx*4 + 3].color = colour;
				}
				else if(sf::Keyboard::isKeyPressed(sf::Keyboard::J))
				{
					double brightness = calculate_curl(x, y);
					brightness = abs(brightness) * 100.0;
					brightness = brightness / (1.0 + brightness);
					sf::Color colour = gradient(brightness, grad);

					pixels[idx*4 + 0].color = colour;
					pixels[idx*4 + 1].color = colour;
					pixels[idx*4 + 2].color = colour;
					pixels[idx*4 + 3].color = colour;
				}
				else
				{
					double u = cells[idx].velocity.x * 5.f;
					double v = cells[idx].velocity.y * 5.f;
					double brightness = sqrt(u*u + v*v);

					//brightness = sqrt(u*u + v*v) / 50.f;
					//brightness = 1.f - (u*u + v*v) / 100.f;

					//float r = brightness * 5.f;
					//float g = brightness * 5.f;
					//float b = brightness * 5.f;

					//u = abs(u);
					//v = abs(v);

					//r += u;
					//g += v;
					//b += sqrt(u*u + v*v);

					//r = r / (1.f + r);
					//g = g / (1.f + g);
					//b = b / (1.f + b);

					//r *= 255.f;
					//g *= 255.f;
					//b *= 255.f;

					brightness = brightness / (1.0 + brightness);
					sf::Color colour = gradient(brightness, grad);

					//sf::Color colour(
					//	(uint8_t)r,
					//	(uint8_t)g,
					//	(uint8_t)b
					//);

					pixels[idx*4 + 0].color = colour;
					pixels[idx*4 + 1].color = colour;
					pixels[idx*4 + 2].color = colour;
					pixels[idx*4 + 3].color = colour;
				}
			}
		}

		memcpy(prev_pixels, pixels, 4 * MAP_WIDTH * MAP_HEIGHT * sizeof(sf::Vertex));

		if (sf::Keyboard::isKeyPressed(sf::Keyboard::S))
			smoothify_blocks(pixels, prev_pixels);

		window.clear();
		window.draw(pixels, 4 * MAP_WIDTH * MAP_HEIGHT, sf::Quads);
		window.display();
	}

	delete[] cells;
	delete[] cells_prev;
	delete[] divergence;
	delete[] p;
	delete[] pixels;
	delete[] prev_pixels;

	return 0;
}
