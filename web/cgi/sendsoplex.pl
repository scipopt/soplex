#!/softis/bin/perl
# $Id: sendsoplex.pl,v 1.2 2002/01/16 14:49:17 bzfkocht Exp $

# einige Vereinbarungen, fuer die Lage von Dateien
# Arbeitsverzeichnis, ggf. aendern
$DOPATH    = "/zibis/Optimization/Software/Soplex/Test/cgi"; 
$COUNTFILE = $DOPATH."/usercount.dat";      # Zaehlerdatei
$DATAFILE  = $DOPATH."/users.dat";          # Benutzerdaten

# in FORMFIELDS stehen die Namen der Eingabefelder
# in FIELDDATA stehen spaeter die zugehoerigen Daten

@FORMFIELDS = 
("tname", "email", "street", "city", "zip", "country", "accept1", "accept2");

# in jedem Datensatz werden die einzelnen Felder durch FSEPARATOR
# getrennt

$FSEPARATOR = ":";

# WWW::CGI ist der Parser, Util::Lock lockt eine Datei

push(@INC,$DOPATH);
require WWW::CGI;
require Util::Lock;
require "pwd.pl";

###################
# Fehlerbehandlung
###################

sub ReadCGI_Fail {
    print "Content-type: text/html\n<html><head>";
    print "<title>CGI-Failure</title></head>";
    print "<body><h3>Server could not read data</h3></body></html>";
    exit(0);
}

sub NoLock_Fail {
    print "Content-type: text/html\n\n<html><head>";
    print "<title>File-Failure</title></head>";
    print "<body><h3>Sorry, somebody else is currently requesting SoPlex. Please try again.</h3></body></html>";
    exit(0);
}

sub CheckFields {
    if ($cgi->{"tname"} eq '') { $fault=1; $FString=$fault.'. Full name<br>'; }
    if (($cgi->{"email"} eq '')|| (!($cgi->{"email"} =~ /.+@.+\..+/))) { 
        $fault++; 
        $FString.=$fault.'. Email<br>'; 
    }
    if ($cgi->{"street"} eq '') { $fault++; $FString.=$fault.'. Street address<br>'; }
    if ($cgi->{"city"} eq '') { $fault++; $FString.=$fault.'. City<br>'; }
    if ($cgi->{"zip"} eq '') { $fault++; $FString.=$fault.'. Zip Code<br>'; }
    if ($cgi->{"country"} eq '') { $fault++; $FString.=$fault.'. Country<br>'; }

    if ($fault) {
        print "Content-type: text/html\n\n<html><head>";
        print "<title>SoPlex Registration</title></head>";
        print "<body><h3>Please fill in the following fields correctly:</h3>".$FString."</body></html>";
    exit(0);
    }
}

sub CheckCertification {
    if ($cgi->{"accept2"} ne 'on') {
        $fault=1;
        $FString="1. You must certify to use the software for noncommercial and academic purposes only.<br>";
    }
    if ($cgi->{"accept1"} ne 'on') {
        $fault++;
        $FString.=$fault.". You must notify and accept the ZIB ACADEMIC LICENSE.";
    }
    if ($fault) {
        print "Content-type: text/html\n\n<html><head>";
        print "<title>SoPlex Registration</title></head>";
        print "<body><h3>Missing Certification:</h3>".$FString."</body></html>";
        exit(0);
    }
}

sub MailSoPlex {
   #system("/usr/ucb/mail -s SoPlex ".@_[0]." < ".$DOPATH."/soplex120.tgz.uu");
    system("sh ".$DOPATH."/sendit.sh ".@_[0]);
}


##########################################################################
# Hier ist der eigentliche Einstiegspunkt des Skripts
#####################################################

$cgi = new WWW::CGI(); # Der Parser
$lock = new Util::Lock($COUNTFILE); # Lockt das Counter-File

# lese
$cgi->read() || ReadCGI_Fail();
$fault=0; # Fehler testen
$FString='';

# Checkboxes und Eingabefelder abfragen
CheckCertification();
CheckFields();

# Locke die Zaehlerdatei
$lock->lock() || NoLock_Fail(); 

# Zaehlerstand einlesen
open (COUNT, $COUNTFILE); 
$counter = <COUNT>;
close (COUNT);

# Zaehlerstand schreiben
open (COUNT, ">".$COUNTFILE);
if (defined($counter)) { $counter++; print COUNT +$counter; }
else { 
    $counter=1;
    print COUNT $counter; 
}
close (COUNT);

###################################################################
# Der Counter konnte inkrementiert werden, schreibe die Userdaten
# (Datensatz einfach anhaengen)
################################

open (DATA, ">>".$DATAFILE);

foreach $elem (@FORMFIELDS) { print DATA $cgi->{$elem},$FSEPARATOR; }

# der Counter wird ans Ende gehaengt

print DATA $counter;
print DATA "\n";
close (DATA);

MailSoPlex($cgi->{"email"});

# Gib den Zaehler wieder frei
$lock->unlock();

###########################################################################
#
# Rueckmeldung an den User
#
##########################

print "Content-type: text/html\n\n";
print "<html><head><title>SoPlex Registration</title></head><body>";
print "Your request was succesful. You are registered with the following data:<br><br>";
print $cgi->{"tname"}."<br>";
print $cgi->{"email"}."<br>";
print $cgi->{"street"}."<br>";
print $cgi->{"city"}."<br>";
print $cgi->{"zip"}."<br>";
print $cgi->{"country"}."<br>";
print "</body></html>";

exit(0);
