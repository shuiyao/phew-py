int InitIons();

double IonFrac(float temp, float density, int ionid);
void load_fraction_tables();
#ifdef NONEQUIL
float Noneq_Frac(float T, float Z, int ionid);
int calc_noneq_fractions(float temp, float Z, int ionid, float *ionfraction);
void load_noneq_fraction_tables();
#endif
