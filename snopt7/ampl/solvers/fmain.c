char **xargv;
extern void MAIN__(void);
main(int argc, char **argv)
{
  xargv = argv;
  MAIN__();
  return 0;
}
