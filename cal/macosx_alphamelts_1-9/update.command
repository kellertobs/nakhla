#!/usr/bin/perl -w

use strict;

my($windows, $pwd, $oldpwd,  $bin, $install, $home, $path, $program, $file, @files);

$windows = '';
if(exists($ENV{'OS'})) {
  if($ENV{'OS'} =~ /Windows/i) {
    $windows = 1;
  }
}

my $done = '';
until ($done) {

    $pwd = $0; # $0 is the name (and path) of the install script
    $pwd =~ s/update.command//;
    $pwd = './' unless ($pwd); # for 'perl ...'
    # for double-clicking on Mac or 'Open with' on Windows
    unless (chdir "$pwd") { 
	warn "Could not change to directory $pwd: $!";
	last;
    }
    
    if($windows) {
	
	$home = "$ENV{'USERPROFILE'}";
	$home =~ s/^\"//;
	$home =~ s/\"$//;

	print "\nUSERPROFILE=$ENV{'USERPROFILE'}\n";
	print "PATH=$ENV{'PATH'}\n\n";

	$bin = "$home\\Documents\\bin";

	$file = "alphamelts.bat run_alphamelts.command.bat";
	# where.exe is not available on XP and seems to hang if STDERR is redirected
	@files = `for \%i in ($file) do \@find /v "echo off" \"\%~\$PATH:i\" 2> nul`;
	
	@files = grep {
	    if (/.+/) {
		chomp;
		s/(^\-+ )(.*)/$2 \-\>/ || s/(\.exe\"*)/$1\r\n/ || s/\%\*/\r\n/;
	    }
	} @files;

    }
    else {

	# $ENV{'PWD'} may be empty, if running as sudo, or be home directory, if double-clicked on Mac
	# $pwd may be './' if not double-clicked.
	$pwd = `pwd`;
	chomp $pwd;
	$home = $ENV{'HOME'};

	print "\nHOME=$ENV{'HOME'}\n";
	print "PATH=$ENV{'PATH'}\n\n";

	$bin = "$home\/bin";

	$file = "alphamelts run_alphamelts.command";
	@files = `which $file 2> /dev/null | xargs ls -go`;
	@files = grep s/(.* \d\d\:\d\d)(.*) \-\> (.*)/$2 \-\> $3/, @files;

    }

    if (@files) {

	print "Currently installed version:\n\n @files\n";
	print "\nContinue with update (y or n)? ";
	$_ = <STDIN>;
	chomp;
	last if (/^$/ || /^n/i);
	print "\n";

    }

    print "Enter full path for installation directory, or press return to use the current directory".
	" given in brackets\n[$pwd]: ";
    $files[0] = <STDIN>;
    print "\nEnter full path of directory to put links in, or press return to use the default location".
	" given in brackets\n[$bin]: ";
    $files[1] = <STDIN>;

    $install = $pwd;
    for (my $i = 0; $i <=1; $i++) {

	$_ = $files[$i];
	unless (/^$/) {
	    
	    s/\s+$//;  # Trim any trailing white space
	    s/\'/\"/g; # Windows can only deal with name enclosed in double quotes;
	    s/^\"//;   # *nix can take single or double
	    s/\"$//;   # Trim any leftover quotes from start and end of list

	    # Fix any apostrophes in the name (if had single quotes or no quotes or double quotes originally)
	    s/\"\\\"\"/\'/g || s/\\\"/\'/g || s/\"/\'/g;	    
	    s/\\//g unless ($windows); # Remove escaping from any other characters (e.g. '(' and ')')

	    ($i) ? $bin = $_ : $install = $_;
	    
	}

    }

    if($windows) {
	$path = `echo ;%PATH%; | find /C /I \";$bin;\"`;
	chomp $path;
    }
    else {
	$path = ($ENV{'PATH'} =~ /^.*:*$bin:*.*/);
    }
    unless (exists $ENV{'SUDO_USER'}) {
	print "\nThe '$bin' directory appears".(($path) ? " " : " not ")."to be in your path.\n".
	    (($path) ? "" : "Please try running install.command instead.\n");
    }

    print "\nContinue with update (y or n)? ";
    $_ = <STDIN>;
    chomp;
    last if (/^$/ || /^n/i);
    print "\n";

    unless ((-d "$install") || (mkdir "$install")) { 
	warn "Could not make directory '$install': $!";
	last;
    }
    unless ((-d "$bin") || (mkdir "$bin")) { 
	warn "Could not make directory '$bin': $!";
	last;
    }
    
    $oldpwd = '';
    $program = '';

    if ($windows) {

	unless ((-f "$bin\\run_alphamelts.command.bat") && (-f "$install\\run_alphamelts.command")) { 
	    warn "Warning: cannot find previously installed version for those directories!";
	}

	if ((-f "alphamelts_win32.exe") && (-f "alphamelts_win32.bin")) {
	    unless (rename "alphamelts_win32.exe", "alphamelts_win32_bak.exe") {
		warn "Cannot rename old executable: $!";
		last;
	    }
	}
	if (-f "alphamelts_win32.bin") {
	    unless (rename "alphamelts_win32.bin", "alphamelts_win32.exe") {
		warn "Cannot rename new executable: $!";
		last;
	    }
	}
	if (-f "alphamelts_win32.exe") {
	    $program = "alphamelts_win32.exe";
	    unless (open(BAT, ">$bin\\alphamelts.bat")) { 
		warn "Could not open alphamelts.bat file: $!";
		last;
	    }
	    print BAT "\@echo off\n";
	    print BAT "\"$install\\alphamelts_win32.exe\"\n";
	    close(BAT);
	}
	else {
	    warn "Cannot find alphamelts executable!";
	    last;
	}

    }
    else {

	unless ((-f "$bin/run_alphamelts.command") && (-f "$install/run_alphamelts.command")) { 
	    warn "Warning: cannot find previously installed version for those directories!";
	}

	$oldpwd = exists $ENV{'OLDPWD'} ? $ENV{'OLDPWD'} : '';

	$program = '';
	if (-f "$oldpwd/alphamelts/build/Debug/alphamelts") {
	    $oldpwd = "$oldpwd/alphamelts/build/Debug";
	    $program = 'alphamelts';
	}
	else {

	    $oldpwd = '';

	    # 64-bit by default (just make sure there is no 64-bit file to get 32-bit one)
	    if ((-f 'alphamelts_wsl.bin') || (-f 'alphamelts_wsl')) {
		$program = 'alphamelts_wsl';
	    }
	    elsif ((-f 'alphamelts_linux64.bin') || (-f 'alphamelts_linux64')) {
		$program = 'alphamelts_linux64';
	    }
	    elsif ((-f 'alphamelts_linux32.bin') || (-f 'alphamelts_linux32')) {
		$program = 'alphamelts_linux32';
	    }
	    elsif ((-f 'alphamelts_macosx64.bin') || (-f 'alphamelts_macosx64')) {
		$program = 'alphamelts_macosx64';
	    }
	    elsif ((-f 'alphamelts_macosx32.bin') || (-f 'alphamelts_macosx32')) {
		$program = 'alphamelts_macosx32';
	    }
	    
	    if ((-f $program) && (-f "$program.bin")) {
		unless (rename $program, "$program\_bak") {
		    warn "Cannot rename old executable: $!";
		    last;
		}
	    }
	    if (-f "$program.bin") {
		unless (rename "$program.bin", $program) {
		    warn "Cannot rename new executable: $!";
		    last;
		}
	    }

	}

	$done = 1; # perl's symlink function has no force
	$file = 'with-readline';
	if (-f "/usr/local/bin/$file") {
	    $done = ((-f "$bin/$file") && !(-l "$bin/$file")) ? '' :
		((system "ln -sf /usr/local/bin/$file $bin/$file") ? '' : $done);
	}
	elsif (-f "/usr/bin/$file") {
	    $done = ((-f "$bin/$file") && !(-l "$bin/$file")) ? '' :
		((system "ln -sf /usr/bin/$file $bin/$file") ? '' : $done);
	}

	# On Mac force will replace a regular file with a link to itself!
	if ($program) {
            $done = ((-f "$bin/alphamelts") && !(-l "$bin/alphamelts")) ? '' : $done;
	    if ($oldpwd && $done) {
		$done = (system "ln -sf $oldpwd/$program $bin/alphamelts") ? '' : $done;
	    }
	    elsif ($done) {
		$done = (system "ln -sf $install/$program $bin/alphamelts") ? '' : $done;
	    }
	    unless ($done) {
		warn "Could not make all links!\n";
		last;
	    }
	}
	else {
	    warn "Cannot find alphamelts executable!";
	    last;
	}
	
    }

    if (!$oldpwd) {
	my $copy = ($windows) ? 'copy /Y' : 'cp -f'; # perl's copy requires a module to be installed
	$done = 1;
	if ($install ne $pwd) {
	    $done = (system "$copy $program \"$install\"") ? '' : $done;
	}
	unless ($windows) {
	    $done = (chmod 0755, "$install/$program") ? $done : '';
	}
	unless ($done) {
	    warn "Could not copy executable file!\n";
	    last;
	}
    }

    print "Update complete!\n\n";

}

print "\nUpdate aborted or incomplete!\n" unless ($done);
print "Press return to finish.\n";
$_ = <STDIN>;
