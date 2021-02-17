//
// test-async-term.cpp
// ~~~~~~~~~
// starts a soplex process and interrupts it from a different thread using the interrupt pointer
// if everything works SoPLEX should abort with [time limit reached].
// the wait time on timer2 might have to be adjusted depending on speed of the machine.
// needs to be compiled with GMP=true BOOST=true BOOST_LIBDIR=<PATH_TO_BOOST_LIBDIR>
//

#include <iostream>
#include <boost/asio.hpp>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "soplex.h"

class printer
{
public:
  printer(boost::asio::io_context& io)
    : timer1_(io, boost::posix_time::seconds(1)),
      timer2_(io, boost::posix_time::seconds(1)),
      interrupt_(false)
  {
    timer1_.async_wait(
          boost::bind(&printer::runSoplex, this));

    timer2_.async_wait(
          boost::bind(&printer::setInterrupt, this));
  }

  void runSoplex()
  {
     soplex::SoPlex soplex;

     soplex.readFile("../check/instances/scagr25.mps");

     soplex.optimize(&interrupt_);

  }
  void doInterrupt()
  {
     interrupt_ = true;
  }

  void setInterrupt()
  {
      timer2_.expires_at(timer2_.expires_at() + boost::posix_time::millisec(15));

      timer2_.async_wait(
            boost::bind(&printer::doInterrupt, this));
  }

private:
  boost::asio::deadline_timer timer1_;
  boost::asio::deadline_timer timer2_;
  bool interrupt_;
};

int main()
{
  boost::asio::io_context io;
  printer p(io);
  boost::thread t(boost::bind(&boost::asio::io_context::run, &io));
  io.run();
  t.join();

  return 0;
}