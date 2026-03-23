# Getting Started with s5cmd:

> **s5cmd is a very fast S3 and local filesystem execution tool (32x faster than s3cmd!).**

* **GitHub:** [**https://github.com/peak/s5cmd**](https://github.com/peak/s5cmd)

## Installation:

* Install `s5cmd` into your conda (miniforge preferrably) base environment on Anemone, Archer2 or JASMIN using:

```bash
# Activate base environment:
source /my/miniforge/bin/activate

# Install from conda-forge channel:
conda install s5cmd
```

## Set-Up:

### 1. Credentials
* To provide your static credentials, you must first create a `credentials` file:

```bash
# Move to root directory:
cd ~
# Create a .aws directory & credentials file:
mkdir ~/.aws
touch credentials

# Update owner permissions:
chmod 600 credentials
```

* In the `credentials` file, we can define multiple `profiles` for each JASMIN object store tenancy we need to access. The name of the profile will be used later when we define an `alias` in our `~/.bashrc`.
* For each profile, we provide our access key (`aws_access_key_id`) and secret key (`aws_secret_access_key`) as follows:

```bash
[noc_msm]
aws_access_key_id="my-noc-msm-access-key"
aws_secret_access_key="my-noc-msm-secret-key"

[rapid_evo]
aws_access_key_id="my-rapid-evolution-access-key"
aws_secret_access_key="my-rapid-evolution-secret-key"
```

### 2. Aliases

* Next we will add two `aliases` to our `~/.bashrc` file to make working with `s5cmd` on the command line more efficient.

* In your `~/.bashrc` file add the following lines...

```bash
alias s3_msm='s5cmd --endpoint-url https://noc-msm-o.s3-ext.jc.rl.ac.uk --profile noc_msm'

alias s3_revo='s5cmd --endpoint-url https://rapidevolution-o.s3-ext.jc.rl.ac.uk --profile rapid_evo'
```

* ...this will map the `s3_msm` command to the base `s5cmd` command which specified the endpoint URL of the `noc-msm` tenancy in the JASMIN object store and provides the necessary credentials to access this tenancy by passing the name of the profile we added to our `~/.aws/credentials` file above.

> **Note:** Using `s3` in the aliases above is a personal preference since we are working with S3-compatible object stores, `s5` would equally be applicable since we are using `s5cmd` rather than `s3cmd` provided by AWS themselves.

## First Usage:

* To ensure you have correctly configured `s5cmd`, we can now list all of the buckets available in the `noc-msm` tenancy of the JASMIN object store as follows:

```bash
s3_msm ls

# -- Returns -- #
2024/11/05 11:41:31  s3://npd12-j001-2005
2024/07/17 15:18:33  s3://mamma-mia
2024/05/31 10:14:03  s3://tobfer6
2023/03/30 13:25:45  s3://senemov2
2021/11/24 13:12:27  s3://NEMO.eORC025
2025/10/08 20:11:44  s3://oceandatastore
2025/06/27 08:46:00  s3://nemotest101
2025/03/03 14:40:15  s3://npd-eorca1-era5v1
2025/03/12 21:10:42  s3://npd-eorca12-era5v1
2024/11/06 13:05:00  s3://tobias-chunk-test20
...
```

* For a complete list of commands available in `s5cmd`, we can use:

```bash
s3_msm

# -- Returns -- #
NAME:
   s5cmd - Blazing fast S3 and local filesystem execution tool

USAGE:
   s5cmd [global options] command [command options] [arguments...]

COMMANDS:
   ls              list buckets and objects
   cp              copy objects
   rm              remove objects
   mv              move/rename objects
   mb              make bucket
   rb              remove bucket
   select          run SQL queries on objects
   du              show object size usage
   cat             print remote object content
   pipe            stream to remote from stdin
   run             run commands in batch
   sync            sync objects
   version         print version
   bucket-version  configure bucket versioning
   presign         print remote object presign url
   head            print remote object metadata
   help, h         Shows a list of commands or help for one command

GLOBAL OPTIONS:
   --credentials-file value       use the specified credentials file instead of the default credentials file
   --dry-run                      fake run; show what commands will be executed without actually executing them (default: false)
   --endpoint-url value           override default S3 host for custom services [$S3_ENDPOINT_URL]
   --help, -h                     show help (default: false)
   --install-completion           get completion installation instructions for your shell (only available for bash, pwsh, and zsh) (default: false)
   --json                         enable JSON formatted output (default: false)
   --log value                    log level: (trace, debug, info, error) (default: info)
   --no-sign-request              do not sign requests: credentials will not be loaded if --no-sign-request is provided (default: false)
   --no-verify-ssl                disable SSL certificate verification (default: false)
   --numworkers value             number of workers execute operation on each object (default: 256)
   --profile value                use the specified profile from the credentials file
   --request-payer value          who pays for request (access requester pays buckets)
   --retry-count value, -r value  number of times that a request will be retried for failures (default: 10)
   --stat                         collect statistics of program execution and display it at the end (default: false)
   --use-list-objects-v1          use ListObjectsV1 API for services that don't support ListObjectsV2 (default: false)
```