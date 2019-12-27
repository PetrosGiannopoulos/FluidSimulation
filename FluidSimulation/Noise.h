#pragma once

#include <random>
#include <cmath>
#include "Vectors.h"
#include <cstdlib>
#include <cstdint>

class Noise {

public:

	const float MATH_PI = 3.14159265359f;
	const float MATH_2PI = MATH_PI * 2.0f;
	const float MATH_HALF_PI = MATH_PI / 2.0f;
	const float MATH_DEG_TO_RAD = MATH_PI / 180.0f;
	const float MATH_RAD_TO_DEG = 180.0f / MATH_PI;
	const float MATH_EPSILON = 1e-6f;

	Vector2 m_gradients[256];
	int  m_permutations[256];
public:

	Noise() {

	}

	static inline float lerp(float a, float b, float v)
	{
		return a * (1 - v) + b * v;
	}

	static inline float Smooth(float v)
	{
		return v * v * (3 - 2 * v);
	}

	static inline float Gradient(const Vector3 &orig, const Vector3 &grad, const Vector3 &p)
	{
		return grad.dot(p - orig);
	}

	static inline float Gradient(const Vector2 &orig, const Vector2 &grad, const Vector2 &p)
	{
		return grad.dot(p - orig);
	}

	template <typename Rnd>
	inline Vector2 RandomGradient2D(Rnd &g)
	{
		
		std::uniform_real_distribution<float> _2pi(0.0f, MATH_PI * 2.0f);
		const float angle = _2pi(g);
		return{ cosf(angle), sinf(angle) };
	}

	Noise(int seed)
	{
		
		std::default_random_engine rnd(seed);
		for (auto &g : m_gradients) {
			g = RandomGradient2D(rnd);
		}

		for (int i = 0; i < 256; i++) {
			int j = std::uniform_int_distribution<int>(0, i)(rnd);
			m_permutations[i] = m_permutations[j];
			m_permutations[j] = i;
		}
	}

	Vector2 get_gradient(int x, int y) const
	{
		int idx =
			m_permutations[x & 255] +
			m_permutations[y & 255];
		return m_gradients[idx & 255];
	}

	void get_gradients(Vector2 *origins, Vector2 *grads, float x, float y) const
	{
		float x0f = floorf(x);
		float y0f = floorf(y);
		int x0 = x0f;
		int y0 = y0f;
		int x1 = x0 + 1;
		int y1 = y0 + 1;

		grads[0] = get_gradient(x0, y0);
		grads[1] = get_gradient(x1, y0);
		grads[2] = get_gradient(x0, y1);
		grads[3] = get_gradient(x1, y1);

		origins[0] = { x0f + 0.0f, y0f + 0.0f };
		origins[1] = { x0f + 1.0f, y0f + 0.0f };
		origins[2] = { x0f + 0.0f, y0f + 1.0f };
		origins[3] = { x0f + 1.0f, y0f + 1.0f };
	}

	float get(float x, float y) const
	{
		Vector2 origins[4];
		Vector2 grads[4];

		get_gradients(origins, grads, x, y);
		float vals[] = {
			Gradient(origins[0], grads[0],{ x, y }),
			Gradient(origins[1], grads[1],{ x, y }),
			Gradient(origins[2], grads[2],{ x, y }),
			Gradient(origins[3], grads[3],{ x, y }),
		};

		float fx = Smooth(x - origins[0].x);
		float vx0 = lerp(vals[0], vals[1], fx);
		float vx1 = lerp(vals[2], vals[3], fx);
		float fy = Smooth(y - origins[0].y);
		return lerp(vx0, vx1, fy);
	}

};