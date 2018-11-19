#ifndef GAMMA_H
#define GAMMA_H

class Gamma {

public:
    Gamma(void) {
        for (size_t i=0; i < 256; i++) {
            _gamma_lut[i] = uint16_t(lrint(_linearize_gamma(i)*65535));
        }
    }
    
    inline uint16_t linearize_gamma(unsigned char x)  const {
        return _gamma_lut[x];
    }

private:

    inline double _linearize_gamma(unsigned char x) const {
        const double C_linear = 0.0404482;
        const double S_linear = 12.9232102;
        const double SRGB_a = 0.055;

        double y = x / 255.0;

        if (y < C_linear) {
            return y / S_linear;
        }
        return pow((y + SRGB_a) / (1.0 + SRGB_a), 2.4);
    }

    uint16_t _gamma_lut[256];    
    
};    
    
#endif // GAMMA_H
