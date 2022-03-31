#!/usr/bin/perl -w

use strict;

my($windows, $pwd, $oldpwd, $bin, $install, $examples, $home, $path, $oldpath, $program, $file, @files, $in_file);

$windows = '';
if(exists($ENV{'OS'})) {
  if($ENV{'OS'} =~ /Windows/i) {
    $windows = 1;
  }
}

my $done = '';
until ($done) {

    $pwd = $0; # $0 is the name (and path) of the install script
    $pwd =~ s/install.command//;
    $pwd = './' unless ($pwd); # for 'perl ...'
    # for double-clicking on Mac or 'Open with' on Windows
    unless (chdir "$pwd") { 
	warn "Could not change to directory $pwd: $!";
	last;
    }
    
    $oldpath = '';
    if($windows) {
	
	$home = "$ENV{'USERPROFILE'}";
	$home =~ s/^\"//;
	$home =~ s/\"$//;

	print "\nUSERPROFILE=$ENV{'USERPROFILE'}\n";
	print "PATH=$ENV{'PATH'}\n\n";

	$bin = "$home\\Documents\\bin";

	$file = "alphamelts.bat run_alphamelts.command.bat column_pick.command.bat file_format.command.bat";
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

	$file = "alphamelts run_alphamelts.command column_pick.command file_format.command";
	@files = `which $file 2> /dev/null | xargs ls -go`;
	@files = grep s/(.* \d\d\:\d\d)(.*) \-\> (.*)/$2 \-\> $3/, @files;

    }

    if (@files) {

	print "Currently installed version:\n\n @files\n";
	print "\nContinue with installation (y or n)? ";
	$_ = <STDIN>;
	chomp;
	last if (/^$/ || /^n/i);
	print "\n";

    }

    $install = $pwd;
    $examples = $pwd.(($windows) ? "\\" : "/")."examples";
    $examples =~ s/\\\\/\\/g;
    $in_file = 'alphamelts_default_env.txt';

    print "Enter full path for installation directory, or press return to use the current directory".
	" given in brackets\n[$install]: ";
    $files[0] = <STDIN>;
    print "\nEnter full path of directory to put links in, or press return to use the default location".
	" given in brackets\n[$bin]: ";
    $files[1] = <STDIN>;
    print "\nEnter full path of directory to put examples in, or press return to use the default location".
	" given in brackets\n[$examples]: ";
    $files[2] = <STDIN>;
    print "\nEnter full or relative path of the default settings file, or press return to use the one".
	" given in brackets\n[$in_file]: ";
    $files[3] = <STDIN>;

    for (my $i = 0; $i <=3; $i++) {

	$_ = $files[$i];
	unless (/^$/) {
	    
	    s/\s+$//;  # Trim any trailing white space
	    s/\'/\"/g; # Windows can only deal with name enclosed in double quotes;
	    s/^\"//;   # *nix can take single or double
	    s/\"$//;   # Trim any leftover quotes from start and end of list

	    # Fix any apostrophes in the name (if had single quotes or no quotes or double quotes originally)
	    s/\"\\\"\"/\'/g || s/\\\"/\'/g || s/\"/\'/g;	    
	    s/\\//g unless ($windows); # Remove escaping from any other characters (e.g. '(' and ')')

	    ($i > 2) ? $in_file = $_ : (($i > 1) ? $examples = $_ : (($i > 0) ? $bin = $_ : $install = $_));
	    
	}
	else {
            $oldpath = 1 if ($i == 1);
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
	    (($path) ? "" : "install.command will try to add it to the path during installation.\n");
    }

    print "\nContinue with installation (y or n)? ";
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
    unless ((-d "$examples") || (mkdir "$examples")) { 
	warn "Could not make directory '$examples': $!";
	last;
    }
    
    $oldpwd = '';
    $program = '';

    if ($windows) {

	# back up relevant parts of registry as .reg file as go along
	# HKLM = HKEY_LOCAL_MACHINE ... settings for all users (requires admin password)
	# HKCU = HKEY_CURRENT_USER ... settings for current user (overrides local machine ones)
	# HKCR = HKEY_CLASSES_ROOT ... a combination of user and machine settings

	my $key = 'HKCU\Software\Classes\.command';
	system "reg export $key command.reg /y" unless (system "reg query $key");
	system "reg add $key /ve /t REG_SZ /d command_auto_file /f";
	print "Command was: reg add $key /ve /t REG_SZ /d command_auto_file /f\n\n";

	$key = 'HKCU\Software\Classes\Applications\perl.exe';
	system "reg export $key perl.reg /y" unless (system "reg query $key");

	# Use whichever perl program has been assigned to open install.command
	system "reg add $key\\shell\\open\\command /ve /t REG_SZ /d \"$^X \\\"%1\\\" %*\" /f";
	print "Command was: reg add $key\\shell\\open\\command /ve /t REG_SZ /d \"$^X \\\"%1\\\" %*\" /f\n\n";

	# Right-clicking allows user to view / edit script
	system "reg add $key\\shell\\edit\\command /ve /t REG_EXPAND_SZ /d".
	    " \"\%SystemRoot\%\\system32\\Notepad.exe %1\" /f";
	print "Command was: reg add $key\\shell\\edit\\command /ve /t REG_EXPAND_SZ /d".
	    " \"\%SystemRoot\%\\system32\\Notepad.exe %1\" /f\n\n";

	$key = 'HKCU\Software\Classes\command_auto_file';
	system "reg export $key command_auto_file.reg /y" unless (system "reg query $key");
	system "reg copy HKCU\\Software\\Classes\\Applications\\perl.exe $key /s /f";
	print "Command was: reg copy HKCU\\Software\\Classes\\Applications\\perl.exe $key /s /f\n\n";


	# 64-bit by default (just make sure there is no 64-bit file to get 32-bit one)
	if ((-f 'alphamelts_win64.bin') || (-f 'alphamelts_win64.exe')) {
	    $program = 'alphamelts_win64';
	}
	elsif ((-f 'alphamelts_win32.bin') || (-f 'alphamelts_win32.exe')) {
	    $program = 'alphamelts_win32';
	}

	if ((-f "$program.exe") && (-f "$program.bin")) {
	    unless (rename "$program.exe", "$program\_bak.exe") {
		warn "Cannot rename old executable: $!";
		last;
	    }
	}
	if (-f "$program.bin") {
	    unless (rename "$program.bin", "$program.exe") {
		warn "Cannot rename new executable: $!";
		last;
	    }
	}

	if ($program) {
	    unless (open(BAT, ">$bin\\alphamelts.bat")) { 
		warn "Could not open alphamelts.bat file: $!";
		last;
	    }
	    print BAT "\@echo off\n";
	    print BAT "\"$install\\$program.exe\"\n";
	    close(BAT);
	}
	else {
	    warn "Warning: installing .command scripts only as cannot find alphamelts executable!";
	}

	@files = ('run_alphamelts', 'column_pick', 'file_format');
	foreach $file (@files) {
	    unless (open(BAT, ">$bin\\$file\.command\.bat")) {
		warn "Could not open $file\.command\.bat file: $!";
		last;
	    }
	    print BAT "\@echo off\n";
	    print BAT "\"$install\\$file\.command\" %*\n";
	    close(BAT);
	}

    }
    else {

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

	@files = ('column_pick.command', 'file_format.command', 'run_alphamelts.command');
	push @files, glob '*.melts';
	push @files, glob '*.txt';
	$file = join ' ', @files;

	if (-x "$pwd/file_format.command") {
	    (exists $ENV{'SUDO_USER'}) ? system "sudo -u $ENV{'SUDO_USER'} $pwd/file_format.command $file" :
		system "$^X $pwd/file_format.command $file";
	}
	else { # previous install means local file_format not executable?
	    (exists $ENV{'SUDO_USER'}) ? system "sudo -u $ENV{'SUDO_USER'} file_format.command $file" :
		system "$^X file_format.command $file";
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
	}
	else {
	    warn "Warning: installing .command scripts only as cannot find alphamelts executable!";
	}

	@files = ('column_pick.command', 'file_format.command', 'run_alphamelts.command');
	for $file (@files) {
	    $done = ((-f "$bin/$file") && !(-l "$bin/$file")) ? '' :
		((system "ln -sf $install/$file $bin/$file") ? '' : $done);
	}

	unless ($done) {
	    warn "Warning: could not make all links! Please check that the links directory is writable.\n";
	    warn "Also, please ensure that the installation and links directories are not the same.\n";
	    last;
	}

    }

    $file = 'run_alphamelts.command';
    unless (rename($file, "$file~")) {
	warn "Could not back up file : $!\n";
	warn "Check that \"$file\" exists and \"$file~\" is not write protected or open elsewhere!\n";
	last;
    }
    unless (open(IN, "<$file~")) {
	warn "Could not open file \"$file~\" : $!\n";
	last;
    }
    unless (open(OUT, ">$file")) {
	warn "Could not open file \"$file\" : $!\n";
	last;
    }
    while (<IN>) {
	s/^\$in_file = .*/\$in_file = '$in_file';/;
	print OUT;
    }
    close(IN);
    close(OUT);
    (unlink "$file~") || (warn "Could not delete backup file \"$file~\" : $!\n");

    my $copy = ($windows) ? 'copy /Y' : 'cp -f'; # perl's copy requires a module to be installed
    @files = glob '*_*.command';
    push @files, $program if ($program && !$oldpwd);

    $done = 1;
    for $file (@files) {
	if ($install ne $pwd) {
	    $done = (system "$copy $file \"$install\"") ? '' : $done;
	}
	unless ($windows) {
	    $done = (chmod 0755, "$install/$file") ? $done : '';
	}
    }
    unless ($done) {
	warn "Could not copy all executable files!\n";
	last;
    }
    
    @files = glob '*.melts *.txt';
    push @files, "column_pick.m";

    $done = 1;
    for $file (@files) {
	$done = (system "$copy $file \"$examples\"") ? '' : $done;
	unless ($windows) {
	    $done = (chmod 0644, "$examples/$file") ? $done : '';
	}
    }
    unless ($done) {
	warn "Could not copy all example files!\n";
	last;
    }

    if ($windows) {

	unless ($path) {

	    my $key = 'HKCU\Environment';

	    $oldpath = ($oldpath) ? `echo ;%PATH%; | find /C /I \";$home\\bin;\"` : '';
            chomp $oldpath;

	    if (system "reg query $key /v Path") {
		print "Adding \"$bin\" to the Path variable.\n";
		system "setx Path \"$bin\"";
		print "Command was: setx Path \"$bin\"\n\n";
	    }
	    else {
		system "reg export $key path.reg /y";
		if ($oldpath) {

		    print "\n        <><><><> IMPORTANT: Updating 'Path' variable! <><><><>\n\n".
                    "As of version 1.5, default folder for links has moved from\n\"$home\\bin\" to \"$bin\".\n\n";
                    print "     <><><><> See http://goo.gl/Wg3n06 for more details. <><><><>\n\n";

		    $oldpath = `reg query $key /v Path | find /i "Path"`;
                    chomp $oldpath;
		    $oldpath =~ s/\s*Path\s+REG_SZ\s+//i;
                    $home =~ s/\\/\\\\/g;
		    $oldpath =~ s/$home\\bin/$home\\Documents\\bin/;
                    $oldpath =~ s/\\\\/\\/g;

		    system "setx Path \"$oldpath\"";
		    print "Command was: setx Path \"$oldpath\"\n\n";
		}
		else {
		    print "\n<><><><> IMPORTANT: the 'Path' variable exists so not going to edit it! <><><><>\n\n".
			"Please ensure that the following directory:\n     $bin\nis in your path...\n\n".
			"On Windows 10, click the Start Button and select 'Settings' or click the cog symbol".
                        "to open Windows Settings. Type 'environment' in the search box and select the option".
			"'Edit environment variables for your account' from the drop-down menu.\n".
			"On Windows 7, find the Control Panel in the Start Menu.\n".
			"Go to: Control Panel->User Accounts and Family Safety->User Accounts.".
			"Select 'Change my environment variables' from the left hand menu.\n\n".
			"Edit the user's Path value. Add the following to the end of the existing data:\n".
			"     ;\"$bin\"\n\n";
                     print "     <><><><> See http://goo.gl/Wg3n06 for more details. <><><><>\n\n";

		    
		}

	    }
	}

    } 
    else {

	unless ($path || exists $ENV{'SUDO_USER'}) {
	    if ((-f "$home/.bash_profile") || (-f "$home/.profile")) {
		my $profile = (-f "$home/.bash_profile") ? "$home/.bash_profile" : "$home/.profile";
		print "\n<><><><> IMPORTANT: the file '$profile' exists so not going to edit it! <><><><>\n\n".
		    "Please ensure that the '$bin' directory is in your path...\n\n".
		    "Edit your '$profile' file (e.g. using 'nano $profile') and include\n".
		    "something like the following:\n\n".
		    "     # set PATH so it includes user's private bin if it exists\n".
		    "     if [ -d \"$bin\" ] ; then\n".
		    "             PATH=\"$bin:\$PATH\"\n".
		    "             export PATH\n".
		    "     fi\n\n".
		    "     <><><><> See http://goo.gl/u0fbfB for more details. <><><><>\n\n";
	    }
	    else {
		unless (open(BASH, ">$home/.profile")) {
		    warn "Could not open .profile file: $!";
		    last;
		}
		print BASH "# ~/.profile: executed by the command interpreter for login shells.\n\n".
		    "# This file is not read by bash(1), if ~/.bash_profile or ~/.bash_login\n".
		    "# exists.\n\n".
		    "# set PATH so it includes user's private bin if it exists\n".
		    "if [ -d \"$bin\" ] ; then\n".
		    "    PATH=\"$bin:\$PATH\"\n".
		    "    export PATH\n".
		    "fi\n";
		print BASH "\nexport DISPLAY=:0\n" if (system "uname -r" =~ /Microsoft/i);
		close(BASH);
	    }
	}

    }

    print "Installation complete!";
    if (exists $ENV{'SUDO_USER'}) {
	print " Please run the install script again without 'sudo'.\n".
	    "Check that the '$bin' directory is in your path." unless ($path);
    }
    else {
	print " Please log out and in and then open this window again.\n".
	    "Check that the '$bin' directory is in your path." unless ($path);
    }
    print "\n\n";

}

print "\nInstallation aborted or incomplete!\n" unless ($done);
print "Press return to finish.\n";
$_ = <STDIN>;
