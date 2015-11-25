#!/usr/bin/perl

#NOTE: This is a utility designed to repair broken SVG files produced by the RHEL-5 cairo package.
# Not intended for general use!
# See https://bugzilla.redhat.com/show_bug.cgi?id=691844 for more information about the exact
# nature of this bug.

#Read line-by-line from stdin or the filename supplied:

while(<>){
  my $line = $_;
  if( $line =~ m/symbol id/ ){
    $line =~ s/symbol id=\"glyph/symbol overflow=\"visible\" id=\"glyph/;
  }
  print $line;
}


