﻿/*
 * Original author: Nicholas Shulman <nicksh .at. u.washington.edu>,
 *                  MacCoss Lab, Department of Genome Sciences, UW
 *
 * Copyright 2011 University of Washington - Seattle, WA
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace pwiz.Common.DataBinding
{
    /// <summary>
    /// Class for exporting data to a file.
    /// </summary>
    public interface IDataFormat
    {
        void WriteRow(TextWriter writer, IEnumerable values);
        string FileFilter { get; }
    }

    public static class DataFormats
    {
        public static readonly IDataFormat Tsv = new TextFormat("\t", "Tab Separated Values(*.tsv)|*.tsv");
        public static readonly IDataFormat Csv = new TextFormat(",", "Comma Separated Values(*.csv)|*.csv");
        class TextFormat : IDataFormat
        {
            public TextFormat(string separator, string fileFilter)
            {
                Separator = separator;
                FileFilter = fileFilter;
            }

            string Separator { get; set; }
            public string FileFilter { get; private set; }
            public void WriteRow(TextWriter writer, IEnumerable values)
            {
                string sep = "";
                foreach (var value in values)
                {
                    writer.Write(sep);
                    sep = Separator;
                    writer.Write(Escape(value));
                }
                writer.WriteLine();
            }
            static readonly Regex RegexQuote = new Regex("\"");
            private string Escape(object o)
            {
                if (o == null)
                {
                    return "";
                }
                return "\"" + RegexQuote.Replace(o.ToString(), "\"\"") + "\"";
            }
        }
    }
}
