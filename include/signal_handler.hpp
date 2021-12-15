#ifndef SIGNAL_HANDLER
#define SIGNAL_HANDLER

#include <signal.h>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>

static bool running = true;

static void ctrl_c_handler (int s)
{
    printf ("Caught signal %d\n", s);
    running = false;
}


static class signal_reg
{

  public:
    signal_reg ()
    {
        struct sigaction sigIntHandler;
        sigIntHandler.sa_handler = ctrl_c_handler;
        sigemptyset (&sigIntHandler.sa_mask);
        sigIntHandler.sa_flags = 0;

        sigaction (SIGINT, &sigIntHandler, NULL);
    }
} _signal_reg;


#endif