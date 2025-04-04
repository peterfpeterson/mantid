import re
import sys
from .cxxtest_misc import abort

# Global variables
suites = []
suite = None
inBlock = 0
options=None
lastLineVoid = False


def scanInputFiles(files, _options):
    '''Scan all input files for test suites'''
    global options
    options=_options
    for file in files:
        scanInputFile(file)
    global suites
    if len(suites) == 0 and not options.root:
        abort( 'No tests defined' )

    return [options,suites]


lineCont_re = re.compile(r'(.*)\\\s*$')


def scanInputFile(fileName):
    '''Scan single input file for test suites'''
    file = open(fileName, encoding='utf-8')
    prev = ""
    lineNo = 0
    contNo = 0
    while 1:
        line = file.readline()
        if not line:
            break
        lineNo += 1

        m = lineCont_re.match(line)
        if m:
            prev += m.group(1) + " "
            contNo += 1
        else:
            scanInputLine( fileName, lineNo - contNo, prev + line )
            contNo = 0
            prev = ""
    if contNo:
        scanInputLine( fileName, lineNo - contNo, prev + line )

    closeSuite()
    file.close()


def scanInputLine( fileName, lineNo, line ):
    '''Scan single input line for interesting stuff'''
    scanLineForExceptionHandling( line )
    scanLineForStandardLibrary( line )

    scanLineForSuiteStart( fileName, lineNo, line )

    global suite
    if suite:
        scanLineInsideSuite( suite, lineNo, line )


def scanLineInsideSuite( suite, lineNo, line ):
    '''Analyze line which is part of a suite'''
    global inBlock
    if lineBelongsToSuite( suite, lineNo, line ):
        scanLineForTest( suite, lineNo, line )
        scanLineForCreate( suite, lineNo, line )
        scanLineForDestroy( suite, lineNo, line )


def lineBelongsToSuite( suite, lineNo, line ):
    '''Returns whether current line is part of the current suite.
    This can be false when we are in a generated suite outside of CXXTEST_CODE() blocks
    If the suite is generated, adds the line to the list of lines'''
    if not suite['generated']:
        return 1

    global inBlock
    if not inBlock:
        inBlock = lineStartsBlock( line )
    if inBlock:
        inBlock = addLineToBlock( suite, lineNo, line )
    return inBlock


std_re = re.compile( r"\b(std\s*::|CXXTEST_STD|using\s+namespace\s+std\b|^\s*\#\s*include\s+<[a-z0-9]+>)" )


def scanLineForStandardLibrary( line ):
    '''Check if current line uses standard library'''
    global options
    if not options.haveStandardLibrary and std_re.search(line):
        if not options.noStandardLibrary:
            options.haveStandardLibrary = 1


exception_re = re.compile( r"\b(throw|try|catch|TSM?_ASSERT_THROWS[A-Z_]*)\b" )


def scanLineForExceptionHandling( line ):
    '''Check if current line uses exception handling'''
    global options
    if not options.haveExceptionHandling and exception_re.search(line):
        if not options.noExceptionHandling:
            options.haveExceptionHandling = 1


classdef = r'(?:::\s*)?(?:\w+\s*::\s*)*\w+'
baseclassdef = r'(?:public|private|protected)\s+%s' % (classdef,)
testsuite = r'(?:(?:::)?\s*CxxTest\s*::\s*)?TestSuite'
suite_re = re.compile( r"\bclass\s+(%s)\s*:(?:\s*%s\s*,)*\s*public\s+%s"
                       % (classdef, baseclassdef, testsuite) )
generatedSuite_re = re.compile( r'\bCXXTEST_SUITE\s*\(\s*(\w*)\s*\)' )


def scanLineForSuiteStart( fileName, lineNo, line ):
    '''Check if current line starts a new test suite'''
    m = suite_re.search( line )
    if m:
        startSuite( m.group(1), fileName, lineNo, 0 )
    m = generatedSuite_re.search( line )
    if m:
        sys.stdout.write( "%s:%s: Warning: Inline test suites are deprecated.\n" % (fileName, lineNo) )
        startSuite( m.group(1), fileName, lineNo, 1 )


def startSuite( name, file, line, generated ):
    '''Start scanning a new suite'''
    global suite
    closeSuite()
    object_name = name.replace(':',"_")
    suite = { 'name'         : name,
              'file'         : file,
              'cfile'        : cstr(file),
              'line'         : line,
              'generated'    : generated,
              'object'       : 'suite_%s' % object_name,
              'dobject'      : 'suiteDescription_%s' % object_name,
              'tlist'        : 'Tests_%s' % object_name,
              'tests'        : [],
              'lines'        : [] }


def lineStartsBlock( line ):
    '''Check if current line starts a new CXXTEST_CODE() block'''
    return re.search( r'\bCXXTEST_CODE\s*\(', line ) is not None


test_re = re.compile( r'^([^/]|/[^/])*\bvoid\s+([Tt]est\w+)\s*\(\s*(void)?\s*\)' )
void_re = re.compile( r'^([^/]|/[^/])*\bvoid\s*$' )
test_novoid_re = re.compile( r'^([^/]|/[^/])*\b([Tt]est\w+)\s*\(\s*(void)?\s*\)' )


def scanLineForTest( suite, lineNo, line ):
    '''Check if current line starts a test'''
    global lastLineVoid
    m = test_re.search( line )
    if m:
        addTest( suite, m.group(2), lineNo )
        lastLineVoid = False
    elif lastLineVoid:
        m2 = test_novoid_re.search(line)
        if m2:
            addTest( suite, m2.group(2), lineNo )
        lastLineVoid = False
    else:
        m3 = void_re.search(line)
        if m3:
            lastLineVoid = True


def addTest( suite, name, line ):
    '''Add a test function to the current suite'''
    test = { 'name'   : name,
             'suite'  : suite,
             'class'  : 'TestDescription_%s_%s' % (suite['object'], name),
             'object' : 'testDescription_%s_%s' % (suite['object'], name),
             'line'   : line,
             }
    suite['tests'].append( test )


def addLineToBlock( suite, lineNo, line ):
    '''Append the line to the current CXXTEST_CODE() block'''
    line = fixBlockLine( suite, lineNo, line )
    line = re.sub( r'^.*\{\{', '', line )

    e = re.search( r'\}\}', line )
    if e:
        line = line[:e.start()]
    suite['lines'].append( line )
    return e is None


def fixBlockLine( suite, lineNo, line):
    '''Change all [E]TS_ macros used in a line to _[E]TS_ macros with the correct file/line'''
    return re.sub( r'\b(E?TSM?_(ASSERT[A-Z_]*|FAIL))\s*\(',
                   r'_\1(%s,%s,' % (suite['cfile'], lineNo),
                   line, 0 )


create_re = re.compile( r'\bstatic\s+\w+\s*\*\s*createSuite\s*\(\s*(void)?\s*\)' )


def scanLineForCreate( suite, lineNo, line ):
    '''Check if current line defines a createSuite() function'''
    if create_re.search( line ):
        addSuiteCreateDestroy( suite, 'create', lineNo )


destroy_re = re.compile( r'\bstatic\s+void\s+destroySuite\s*\(\s*\w+\s*\*\s*\w*\s*\)' )


def scanLineForDestroy( suite, lineNo, line ):
    '''Check if current line defines a destroySuite() function'''
    if destroy_re.search( line ):
        addSuiteCreateDestroy( suite, 'destroy', lineNo )


def cstr( s ):
    '''Convert a string to its C representation'''
    return '"' + s.replace( '\\', '\\\\' ) + '"'


def addSuiteCreateDestroy( suite, which, line ):
    '''Add createSuite()/destroySuite() to current suite'''
    if which in suite:
        abort( '%s:%s: %sSuite() already declared' % ( suite['file'], str(line), which ) )
    suite[which] = line


def closeSuite():
    '''Close current suite and add it to the list if valid'''
    global suite
    if suite is not None:
        if len(suite['tests']) != 0:
            verifySuite(suite)
            rememberSuite(suite)
        suite = None


def verifySuite(suite):
    '''Verify current suite is legal'''
    if 'create' in suite and 'destroy' not in suite:
        abort( '%s:%s: Suite %s has createSuite() but no destroySuite()' %
               (suite['file'], suite['create'], suite['name']) )
    elif 'destroy' in suite and 'create' not in suite:
        abort( '%s:%s: Suite %s has destroySuite() but no createSuite()' %
               (suite['file'], suite['destroy'], suite['name']) )


def rememberSuite(suite):
    '''Add current suite to list'''
    global suites
    suites.append( suite )

#
# Copyright 2008 Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
# retains certain rights in this software.
#
