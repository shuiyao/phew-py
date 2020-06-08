void get_snap_string(int snapnum, char *snapstr);
void read_header(struct dump *header, FILE *ftipsy);
int headerinfo(struct dump *header);
#ifdef IONS
int InitIons();
double IonFrac(float temp, float density, int ionid);
void load_fraction_tables();
#endif
